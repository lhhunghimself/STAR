#ifndef CODE_BAMTagBuffer
#define CODE_BAMTagBuffer

#include "IncludeDefine.h"
#include "BAMoutput.h"
#include "SoloCommon.h"
#include <string>
#include <vector>
#include <mutex>
#include <cstdint>

struct BAMTagEntry {
    uint64_t recordIndex;
    uint32_t readId;    // replaces iReadAll, stores iReadAll cast to 32-bit
    uint16_t alignIdx;
    uint8_t  mate;
    uint8_t  padding;   // for alignment, set to 0
    
    BAMTagEntry() : recordIndex(0), readId(0), alignIdx(0), mate(0), padding(0) {}
    BAMTagEntry(const BAMRecordMeta& meta) 
        : recordIndex(meta.recordIndex), readId(static_cast<uint32_t>(meta.iReadAll))
        , alignIdx(static_cast<uint16_t>(meta.alignIdx)), mate(static_cast<uint8_t>(meta.mate)), padding(0) {}
};

class BAMTagBuffer {
public:
    BAMTagBuffer();
    ~BAMTagBuffer();
    
    // Thread-safe append of alignment metadata
    void append(const BAMRecordMeta& meta);
    
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
    std::vector<BAMTagEntry> entries;
    std::mutex entriesMutex; // Thread safety for append operations
};

#endif
