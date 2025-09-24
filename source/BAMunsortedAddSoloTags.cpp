#include "BAMunsortedAddSoloTags.h"
#include "ErrorWarning.h"
#include "BAMfunctions.h"
#include <fstream>

void BAMunsortedAddSoloTags(const std::string &tmpPath,
                            const std::string &outBamPath,
                            Parameters &P,
                            Genome &genome,
                            Solo &solo)
{
    // Guard: only proceed if CB/UB tags are requested
    if (!solo.pSolo.samAttrYes) {
        // Nothing to inject, just return
        return;
    }

    // Open the tmp file for reading (binary)
    std::ifstream tmpStream(tmpPath.c_str(), std::ios::binary);
    if (!tmpStream.is_open()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open solo tmp file for reading: " << tmpPath << "\n";
        errOut << "SOLUTION: check that the file exists and has proper permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Open the final unsorted BAM for writing
    BGZF *bgzfOut = bgzf_open(outBamPath.c_str(), ("w"+to_string((long long) P.outBAMcompression)).c_str());
    if (bgzfOut == NULL) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not create final unsorted BAM file: " << outBamPath << "\n";
        errOut << "SOLUTION: check that the directory exists and has proper permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Write BAM header immediately
    outBAMwriteHeader(bgzfOut, P.samHeader, genome.chrNameAll, genome.chrLengthAll);

    // Allocate scratch buffers
    char *bam0 = new char[BAMoutput_oneAlignMaxBytes + sizeof(uint64)]; // space for record + trailer
    char *bam1 = new char[BAM_ATTR_MaxSize]; // workspace for addBAMtags

    // Loop until EOF: read records from tmp file and process them
    uint32 size_with_len;
    uint64 trailer;
    
    while (tmpStream.read(reinterpret_cast<char*>(&size_with_len), sizeof(uint32))) {
        // Guardrail: Check record size before reading
        if (size_with_len > BAMoutput_oneAlignMaxBytes) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: BAM record size " << size_with_len << " exceeds maximum buffer size " << BAMoutput_oneAlignMaxBytes << "\n";
            errOut << "SOLUTION: increase BAMoutput_oneAlignMaxBytes or check for corrupted tmp file: " << tmpPath << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        // Read the payload (BAM record)
        if (!tmpStream.read(bam0, size_with_len)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not read BAM record payload from solo tmp file: " << tmpPath << "\n";
            errOut << "SOLUTION: the tmp file may be corrupted\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Read the trailer (contains iReadAll)
        if (!tmpStream.read(reinterpret_cast<char*>(&trailer), sizeof(uint64))) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: could not read trailer from solo tmp file: " << tmpPath << "\n";
            errOut << "SOLUTION: the tmp file may be corrupted\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Append the trailer to the BAM record (addBAMtags expects it there)
        memcpy(bam0 + size_with_len, &trailer, sizeof(uint64));
        
        uint32 size0 = size_with_len; // Initial size (will be updated by addBAMtags)
        
        // Call addBAMtags to inject CB/UB tags
        char *bam0_ptr = bam0; // addBAMtags takes char*& so we need a pointer to modify
        solo.soloFeat[solo.pSolo.featureInd[solo.pSolo.samAttrFeature]]->addBAMtags(bam0_ptr, size0, bam1);

        // Guardrail: Check final size after tag injection
        if (size0 > BAM_ATTR_MaxSize) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: BAM record size after tag injection " << size0 << " exceeds scratch buffer size " << BAM_ATTR_MaxSize << "\n";
            errOut << "SOLUTION: increase BAM_ATTR_MaxSize or check for excessive tag data\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        // Write the resulting record to the final BAM (without trailer)
        bgzf_write(bgzfOut, bam0_ptr, size0);
    }

    // Check if we exited due to EOF or error
    if (!tmpStream.eof()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: error reading from solo tmp file: " << tmpPath << "\n";
        errOut << "SOLUTION: the tmp file may be corrupted\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Cleanup: close files and free memory
    delete[] bam0;
    delete[] bam1;
    
    tmpStream.close();
    
    if (bgzf_close(bgzfOut) != 0) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not close final unsorted BAM file: " << outBamPath << "\n";
        errOut << "SOLUTION: check disk space and permissions\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    // Delete the tmp file on success
    if (remove(tmpPath.c_str()) != 0) {
        // Non-fatal warning - tmp file cleanup failed but processing succeeded
        P.inOut->logMain << "WARNING: could not delete temporary file: " << tmpPath << std::endl;
    }
}
