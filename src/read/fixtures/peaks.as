table peaks
"Peaks as a peak caller writes them: six columns of BED and four of its own."
(
    string chrom;       "The sequence"
    uint   chromStart;  "Where the peak starts, from 0"
    uint   chromEnd;    "Where it ends, the base after its last"
    string name;        "The peak's name"
    uint   score;       "A score from 0 to 1000"
    char[1] strand;     "+, - or ."
    float  signalValue; "How enriched the peak is"
    float  pValue;      "-log10 of its p-value"
    float  qValue;      "-log10 of its q-value"
    int    peak;        "Where its summit is, from its start"
)
