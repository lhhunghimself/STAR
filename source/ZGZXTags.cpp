#include "ZGZXTags.h"
#include <sstream>
#include <algorithm>

string ZGZXTags::formatZGTag(const ReadAnnotFeature &annFeat, const Transcriptome &transcriptome) {
    if (annFeat.fSet.empty()) {
        return "-";
    }
    
    // Convert gene indices to gene IDs and join with commas
    string result;
    bool first = true;
    for (uint32 geneIdx : annFeat.fSet) {
        if (!first) {
            result += ",";
        }
        if (geneIdx < transcriptome.geID.size()) {
            result += transcriptome.geID[geneIdx];
        } else {
            result += "UNKNOWN_GENE_" + to_string(geneIdx);
        }
        first = false;
    }
    
    return result;
}

string ZGZXTags::formatZXTag(const ReadAnnotFeature &annFeat) {
    return overlapTypeToString(annFeat.ovType);
}

const char* ZGZXTags::overlapTypeToString(uint32 ovType) {
    switch (ovType) {
        case ReadAnnotFeature::overlapTypes::none:
            return "none";
        case ReadAnnotFeature::overlapTypes::exonic:
            return "exonic";
        case ReadAnnotFeature::overlapTypes::exonicAS:
            return "exonic";  // Treat antisense exonic as exonic for simplicity
        case ReadAnnotFeature::overlapTypes::exonic50p:
            return "exonic";  // Treat 50% exonic as exonic
        case ReadAnnotFeature::overlapTypes::exonic50pAS:
            return "exonic";  // Treat antisense 50% exonic as exonic
        case ReadAnnotFeature::overlapTypes::intronic:
            return "intronic";
        case ReadAnnotFeature::overlapTypes::intronicAS:
            return "intronic";  // Treat antisense intronic as intronic
        case ReadAnnotFeature::overlapTypes::intergenic:
            return "none";  // Treat intergenic as none
        default:
            return "spanning";  // Default for unknown types
    }
}
