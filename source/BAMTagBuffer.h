#ifndef CODE_BAMTagBuffer
#define CODE_BAMTagBuffer

#include "IncludeDefine.h"
#include "BAMoutput.h"
#include "SoloCommon.h"
#include <string>
#include <vector>
#include <mutex>

struct BAMTagEntry {
    uint64_t recordIndex;
    uint64_t iReadAll;
    uint32_t mate;
    uint32_t alignIdx;
    uint32_t groupIndex;     // Index into readGroups array
    
    BAMTagEntry() : recordIndex(0), iReadAll(0), mate(0), alignIdx(0), groupIndex(UINT32_MAX) {}
    BAMTagEntry(const BAMRecordMeta& meta) 
        : recordIndex(meta.recordIndex), iReadAll(meta.iReadAll), mate(meta.mate), alignIdx(meta.alignIdx)
        , groupIndex(UINT32_MAX) {}
};

class BAMTagBuffer {
public:
    BAMTagBuffer();
    ~BAMTagBuffer();
    
    // Thread-safe append of alignment metadata
    void append(const BAMRecordMeta& meta);
    
    // Hint expected number of reads to pre-size buckets
    void reserveReadCapacity(uint64_t nReads);
    
    // Write tag table to file (called at end of processing)
    void writeTagTable(const std::string& outputPath,
                       const std::vector<readInfoStruct>& readInfo,
                       const std::vector<std::string>& cbWLstr,
                       size_t umiLength);
    
    // Clear buffer to free memory
    void clear();
    
    // Get number of entries
    size_t size() const { return entries.size(); }

private:
    // ReadGroup structure for minimal shared data across alignments of the same read
    struct ReadGroup {
        uint32_t firstEntryIdx;  // index into entries, UINT32_MAX when unused
        uint32_t entryCount;     // number of alignments for this read
        
        ReadGroup() : firstEntryIdx(UINT32_MAX), entryCount(0) {}
    };
    
    uint32_t getOrCreateGroup(uint64_t iReadAll);

    std::vector<BAMTagEntry> entries;
    std::vector<ReadGroup> readGroups;
    std::mutex entriesMutex; // Thread safety for append operations
};

#endif
