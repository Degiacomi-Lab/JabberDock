#!/usr/bin/awk -f

# trim functions for whitespace recognition (in case atomnames are not justified correctly)
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

BEGIN{
    # Reading Fixed-Width Data (see: https://goo.gl/SmjwUt)
    # might need to change last fieldlength (should be 4 but made 6 for extra two spaces)
    FIELDWIDTHS = "6 5 1 4 1 4 1 4 1 3 8 8 8 6 6 6 6"
    # $2: Atom serial number
    # $4: Atom type
    # $5: altLoc; alternate location indicator.
    # $6: Resname
    # $8: ChainID
    # $9: Resid
    # $12: x
    # $13: y
    # $14: z
    # this cuts off the remaining fields (e.g. beta factor)
    #var="   N"
    # ARGV[1] is the script name
    atom_from=trim(ARGV[2])
    atom_to=ARGV[3]
    res_from=ARGV[4]
    res_to=ARGV[5]
    ARGV[2]="" # nuke the variables so it doesn't seek the files
    ARGV[3]=""
    ARGV[4]=""
    ARGV[5]=""
}

{
    #if ($6 == "PE  " && $4 == " N31"){
    #    printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,var,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)
    #}
    if ($6 == res_from && trim($4) == atom_from){
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,atom_to,$5,res_to,$7,$8,$9,$10,$11,$12,$13,$14)
    }
    else {
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)
    }
} 