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
    readGroups.reserve(100000); // Start with 100K read groups capacity
}

BAMTagBuffer::~BAMTagBuffer() {
    clear();
}

uint32_t BAMTagBuffer::getOrCreateGroup(uint64_t iReadAll) {
    // Ensure we have enough space in readGroups
    if (iReadAll >= readGroups.size()) {
        size_t newSize = readGroups.size();
        if (newSize == 0) {
            newSize = 1024;
        }
        
        const uint64_t required = iReadAll + 1;
        while (static_cast<uint64_t>(newSize) < required) {
            if (newSize >= std::numeric_limits<size_t>::max() / 2) {
                newSize = static_cast<size_t>(required);
                break;
            }
            size_t doubled = newSize * 2;
            if (doubled <= newSize) {
                newSize = static_cast<size_t>(required);
                break;
            }
            newSize = doubled;
        }
        
        readGroups.resize(newSize);
    }
    
    uint32_t groupIndex = static_cast<uint32_t>(iReadAll);
    return groupIndex;
}

void BAMTagBuffer::append(const BAMRecordMeta& meta) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    size_t entryIdx = entries.size();
    
    // Get/create read group
    uint32_t groupIndex = getOrCreateGroup(meta.iReadAll);
    
    // Create entry with reference to read group
    BAMTagEntry entry(meta);
    entry.groupIndex = groupIndex;
    
    entries.push_back(entry);
    
    // Update read group bookkeeping
    ReadGroup& group = readGroups[groupIndex];
    if (group.firstEntryIdx == UINT32_MAX) {
        group.firstEntryIdx = static_cast<uint32_t>(entryIdx);
    }
    group.entryCount++;
}


void BAMTagBuffer::reserveReadCapacity(uint64_t nReads) {
    std::lock_guard<std::mutex> lock(entriesMutex);
    if (nReads > readGroups.size()) {
        readGroups.resize(nReads);
    }
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
        if (entry.iReadAll >= readInfo.size()) {
            continue; // Skip invalid entries
        }
        
        const readInfoStruct& readData = readInfo[entry.iReadAll];
        
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
                << entry.iReadAll << "\t"
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
    readGroups.clear();
    readGroups.shrink_to_fit();
}
