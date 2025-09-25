#include "BAMTagBuffer.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <limits>

BAMTagBuffer::BAMTagBuffer() {
    // Reserve some initial capacity to reduce reallocations
    entries.reserve(1000000); // Start with 1M entries capacity
}

BAMTagBuffer::~BAMTagBuffer() {
    clear();
}

void BAMTagBuffer::append(const BAMRecordMeta& meta) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    
    // Guard against overflow of readId field
    if (meta.iReadAll > UINT32_MAX) {
        std::cerr << "ERROR: Read ID " << meta.iReadAll << " exceeds UINT32_MAX limit. "
                  << "Dataset too large for current implementation." << std::endl;
        exit(1);
    }
    
    // Guard against overflow of alignIdx field
    if (meta.alignIdx > UINT16_MAX) {
        std::cerr << "ERROR: Alignment index " << meta.alignIdx << " exceeds UINT16_MAX limit. "
                  << "Too many alignments per read." << std::endl;
        exit(1);
    }
    
    // Guard against overflow of mate field
    if (meta.mate > UINT8_MAX) {
        std::cerr << "ERROR: Mate value " << meta.mate << " exceeds UINT8_MAX limit." << std::endl;
        exit(1);
    }
    
    // Create compact entry
    BAMTagEntry entry;
    entry.recordIndex = meta.recordIndex;
    entry.readId = static_cast<uint32_t>(meta.iReadAll);
    entry.alignIdx = static_cast<uint16_t>(meta.alignIdx);
    entry.mate = static_cast<uint8_t>(meta.mate);
    entry.padding = 0;
    
    entries.push_back(entry);
}



void BAMTagBuffer::writeTagTable(const std::string& outputPath,
                                  const std::vector<readInfoStruct>& readInfo,
                                  const std::vector<std::string>& cbWLstr,
                                  size_t umiLength) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "ERROR: Cannot open tag table output file: " << outputPath << std::endl;
        return;
    }
    
    // Write header without QNAME column
    outFile << "# bam_record_index\tiReadAll\tmate\talign_idx\tCB\tUB\tstatus\n";
    
    // Sort entries by recordIndex to ensure proper ordering
    std::sort(entries.begin(), entries.end(), 
              [](const BAMTagEntry& a, const BAMTagEntry& b) {
                  return a.recordIndex < b.recordIndex;
              });
    
    // Write entries - derive CB/UB/status from readInfo
    for (const auto& entry : entries) {
        if (entry.readId >= readInfo.size()) {
            continue; // Skip invalid entries
        }
        
        const readInfoStruct& readData = readInfo[entry.readId];
        
        // Derive CB from readInfo
        std::string cb = "-";
        if (readData.cb != (uint64)-1 && readData.cb < cbWLstr.size()) {
            cb = cbWLstr[readData.cb];
        }
        
        // Derive UMI from readInfo
        std::string ub = "-";
        if (readData.umi != (uint32)-1) {
            ub = convertNuclInt64toString(readData.umi, umiLength);
        }
        
        // Compute status inline based on presence/absence of CB/UB
        const char* status;
        bool hasCB = (readData.cb != (uint64)-1);
        bool hasUMI = (readData.umi != (uint32)-1);
        
        if (hasCB && hasUMI) {
            status = "OK";
        } else if (!hasCB && !hasUMI) {
            status = "NO_CB_UMI";
        } else if (!hasCB) {
            status = "NO_CB";
        } else {
            status = "NO_UMI";
        }
        
        outFile << entry.recordIndex << "\t"
                << entry.readId << "\t"
                << entry.mate << "\t"
                << entry.alignIdx << "\t"
                << cb << "\t"
                << ub << "\t"
                << status << "\n";
    }
    
    outFile.close();
    
    std::cout << "Tag table written to " << outputPath << " with " << entries.size() << " records (QNAME omitted)" << std::endl;
}

void BAMTagBuffer::clear() {
    std::lock_guard<std::mutex> lock(entriesMutex);
    entries.clear();
    entries.shrink_to_fit();
}
