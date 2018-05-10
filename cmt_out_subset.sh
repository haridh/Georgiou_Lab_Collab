for f in *.txt; do
    pfx=${f%%.txt}; echo "$pfx";
    grep -o '^[^#]*' "$pfx".txt > "$pfx"_cmt_out.txt
    grep "FCGR1" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "FCGR2" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "FCGR3" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "FCRL" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "KIR3DL3" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "LILRA4" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "KIR2DS2" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
    grep "KIR2DS4" "$pfx"_cmt_out.txt | cut -f11,12,19 >> "$pfx"_subset.txt
done
