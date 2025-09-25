#include "SoloFeature.h"
#include "BAMTagBuffer.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"

void SoloFeature::writeTagTableIfRequested(bool filteredPass)
{
    if (!pSolo.writeTagTableEnabled || featureType != pSolo.samAttrFeature) {
        return;
    }
    
    if (!pSolo.bamTagBuffer) {
        P.inOut->logMain << "WARNING: Tag table export requested but BAMTagBuffer is not available" << endl;
        return;
    }
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Writing tag table from BAMTagBuffer to " << pSolo.writeTagTablePath << " (QNAME omitted)" << endl;
    
    // Use BAMTagBuffer to write the tag table with readInfo data
    pSolo.bamTagBuffer->writeTagTable(pSolo.writeTagTablePath, readInfo, pSolo.cbWLstr, pSolo.umiL);
    
    // Clear the buffer to free memory
    pSolo.bamTagBuffer->clear();
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished writing tag table and cleared buffer" << endl;
}