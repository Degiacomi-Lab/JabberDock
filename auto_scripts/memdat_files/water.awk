#!/usr/bin/awk -f

# trim functions for whitespace recognition (in case atomnames are not justified correctly)
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

BEGIN{
    # Reading Fixed-Width Data (see: https://goo.gl/SmjwUt)
    # might need to change last fieldlength (should be 4 but made 6 for extra two spaces)
    FIELDWIDTHS = "6 5 1 4 1 4 1 4 1 3 8 8 8 6 6 6 6"
    O_name="OW  "
    res = "SOL "
    H1_name="HW1 "
    H2_name="HW2 "
}
{
    if ( trim($6)=="WAT" && trim($4)=="O" ) {
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,O_name,$5,res,$7,$8,$9,$10,$11,$12,$13,$14)
    }
    else if ( trim($6)=="WAT" && trim($4)=="H1" ) {
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,H1_name,$5,res,$7,$8,$9,$10,$11,$12,$13,$14)
    }
    else if ( trim($6)=="WAT" && trim($4)=="H2" ) {
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,H2_name,$5,res,$7,$8,$9,$10,$11,$12,$13,$14)
    }
    else {
        printf("%-6s%5s%1s%4s%1s%4s%1s%4s%1s%3s%8s%8s%8s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)
    }
} 