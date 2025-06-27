#!/bin/bash
# maketable.sh: convert gc.prt into const char predefined_tables[N][2][65]
# Usage: ./maketable.sh < gc.prt

# Hardcoded data for private shenanigans
declare -a GAPS_DES GAPS_EAA GAPS_SEA

add_gap() {
  GAPS_DES[$1]="$2"
  GAPS_EAA[$1]="$3"
  GAPS_SEA[$1]="$4"
}

add_gap 0 "11 + Pyl (Kivenson et al. 2021)" \
  "FFLLSSSSYY*OCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" \
  "---M------*---*----M---------------M----------------------------"

add_gap 7 "Unused - DO NOT USE" \
  "***M***************M***************M****************************" \
  "***M***************M***************M****************************"

add_gap 8 "Unused - DO NOT USE" \
  "***M***************M***************M****************************" \
  "***M***************M***************M****************************"

add_gap 17 "Enterosoma/UBA4682 (Shulgina & Eddy 2021)" \
  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRMVVVVAAAADDEEGGGG" \
  "---M------**--*----M---------------M----------------------------"

add_gap 18 "Peptacetobacter (Shulgina & Eddy 2021)" \
  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRQIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" \
  "---M------**--*----M---------------M----------------------------"

add_gap 19 "Anaerococcus, Onthovivens/UBA4855 (Shulgina & Eddy 2021)" \
  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRWIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" \
  "---M------**--*----M---------------M----------------------------"

add_gap 20 "Absconditabacterales (Shulgina & Eddy 2021)" \
  "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRWWIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" \
  "---M------**-----------------------M---------------M------------"

declare -a NCBI_DES NCBI_EAA NCBI_SEA
maxid=0

SYMBOLS=("::=" "{" "}" "," "name" "id" "ncbieaa" "sncbieaa" "Genetic-code-table")
for sym in "${SYMBOLS[@]}"; do
  SYMBOLS_REGEX+="|"
  if [[ $sym == "{" || $sym == "}" ]]; then
    SYMBOLS_REGEX+="\\$sym"
  else
    SYMBOLS_REGEX+="$sym"
  fi
done
SYMBOLS_REGEX="(${SYMBOLS_REGEX:1})"

shopt -s extglob
TOKENS=()
while read -r line; do
  linelen=${#line}
  ((i = 0))
  while ((i < linelen)); do
    slice="${line:i}"
    if [[ $slice =~ ^[$' \t']+ ]]; then
      ((i += ${#BASH_REMATCH[0]}))
      slice="${line:i}"
    fi
    if [[ $slice =~ ^$SYMBOLS_REGEX ]]; then
      TOKENS+=("${BASH_REMATCH[0]}")
      ((i += ${#BASH_REMATCH[0]}))
    elif [[ $slice =~ ^[0-9]+ ]]; then
      TOKENS+=($'\t'"${BASH_REMATCH[0]}")
      ((i += ${#BASH_REMATCH[0]}))
    elif [[ $slice =~ ^\"([^\"]+)\" ]]; then
      # mark string with \n
      TOKENS+=($'\n'"${BASH_REMATCH[1]}")
      ((i += ${#BASH_REMATCH[0]}))
    elif [[ $slice =~ ^\"([^\"]+)$ ]]; then
      # incomplete string
      if ! read -r line2; then
        printf 'Unexpected end of file, want closing \"' >&2
        exit 1
      fi
      line+=" $line2"
      linelen=${#line}
      continue
    elif [[ $slice =~ ^-- ]]; then
      continue 2
    else
      # shellcheck disable=2319
      ret=$?
      printf "Failed to tokenize (%d): %q\n" $ret "$slice" >&2
      exit 1
    fi
  done
done

STATE=FILE_BEGIN
NTOK=${#TOKENS[@]}
((itok = 0))

consume_tok() {
  tok="${TOKENS[$itok]}"
  ((itok++))
  if ((itok >= NTOK)); then
    tok="EOF"
  fi
}
die_unexpected() {
  local state=$STATE
  if [[ -n $1 ]]; then
    state+=" ($1)"
  fi
  printf "Unexpected token %q (%d) in table given %s\n" "$tok" "$itok" "$state" >&2
  exit 1
}
want_tok() {
  # shellcheck disable=2053
  if [[ $tok != $1 ]]; then
    die_unexpected "want $1"
  fi
}
want_string() {
  if [[ $tok != $'\n'* ]]; then
    die_unexpected "want string"
  fi
  tok="${tok:1}"
}
want_num() {
  if [[ $tok != $'\t'* ]]; then
    die_unexpected "want number"
  fi
  tok="${tok:1}"
}

for ((itok = 0; itok < NTOK; )); do
  consume_tok
  case $STATE in
    FILE_BEGIN)
      want_tok "Genetic-code-table"
      consume_tok
      want_tok "::="
      consume_tok
      want_tok "{"
      STATE=TABLE_EXPECT_ELEMENT
      ;;
    TABLE_EXPECT_ELEMENT)
      want_tok "{"
      STATE=ELEMENT_EXPECT_NAME
      ;;
    ELEMENT_EXPECT_NAME)
      want_tok "name"
      consume_tok
      want_string
      name=$tok
      consume_tok
      want_tok ,

      STATE=ELEMENT_EXPECT_ID
      consume_tok
      while [[ $tok == "name" ]]; do
        consume_tok
        want_string
        consume_tok
        want_tok ,
        consume_tok
      done
      want_tok id
      consume_tok
      want_num
      id=$tok
      ((maxid = id > maxid ? id : maxid))
      consume_tok
      want_tok ,

      STATE=ELEMENT_EXPECT_EAA
      consume_tok
      want_tok "ncbieaa"
      consume_tok
      want_string
      eaa=$tok
      consume_tok
      want_tok ,

      STATE=ELEMENT_EXPECT_SEA
      consume_tok
      want_tok "sncbieaa"
      consume_tok
      want_string
      sea=$tok

      STATE=ELEMENT_EXPECT_END
      consume_tok
      want_tok "}"
      NCBI_DES[id]=$name
      NCBI_EAA[id]=$eaa
      NCBI_SEA[id]=$sea

      STATE=ELEMENT_END
      consume_tok
      case $tok in
        EOF) STATE=TABLE_END ;;
        ,) STATE=TABLE_EXPECT_ELEMENT ;;
        *) die_unexpected "want , or }" ;;
      esac
      ;;
    TABLE_END)
      want_tok EOF
      ;;
    *)
      die_unexpected "unknown state"
      ;;
  esac
done

# Make sure initialization only happens at ATG, GTG, TTG
SEAREF="---M---------------M---------------M----------------------------"
mask_sea() {
  local oldsea=$1 newsea
  for i in {0..63}; do
    if [[ ${oldsea:i:1} != '*' && ${oldsea:i:1} != '-' ]]; then
      newsea+=${SEAREF:i:1}
    else
      newsea+=${oldsea:i:1}
    fi
  done
  echo "$newsea"
}

exec >table-data.hh
cat <<EOF
/* This is not a C++ file, but we don't want make glob to see it! */
/* Generated by maketable.sh */
#ifndef TABLE_DATA_HH
#define TABLE_DATA_HH
#define MAXTABLE $maxid
#define STRMAXTABLE "$maxid"
#endif
EOF

exec >table-data.cc
cat <<EOF
/* This is not a C++ file, but we don't want make glob to see it! */
/* Generated by maketable.sh */
#ifndef TABLE_DATA_CC
#define TABLE_DATA_CC
const char predefined_tables[$((maxid + 1))][2][65] = {
EOF

for ((i = 0; i <= maxid; i++)); do
  if [[ -n ${NCBI_DES[i]} ]]; then
    if [[ -n ${GAPS_DES[i]} ]]; then
      printf 'Warning: Both NCBI and GAPS entries for id %d, using NCBI\n' "$i" >&2
    fi
  else
    if [[ -n ${GAPS_DES[i]} ]]; then
      NCBI_DES[i]="(NON-NCBI) "${GAPS_DES[i]}
      NCBI_EAA[i]=${GAPS_EAA[i]}
      NCBI_SEA[i]=${GAPS_SEA[i]}
    else
      NCBI_DES[i]="(NON-NCBI) Undefined, fill with standard"
      NCBI_EAA[i]="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
      NCBI_SEA[i]="---M------**--*----M---------------M----------------------------"
    fi
  fi
  printf '  /* %d. %s */\n  {\n    "%s",\n    "%s"\n  },\n' "$i" "${NCBI_DES[$i]}" "${NCBI_EAA[$i]}" "$(mask_sea "${NCBI_SEA[$i]}")"
done

printf '};\n#endif\n'
