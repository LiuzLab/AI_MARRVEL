process PARSE_TRIO_PED {
    input:
    path ped

    output:
    path "sample_ids.env", emit: sample_ids

    script:
    """
    if [ ! -s "$ped" ]; then
        echo "ERROR: PED file '$ped' is empty or missing." >&2
        exit 1
    fi

    row_count=\$(awk 'NF' "$ped" | wc -l)
    if [ "\$row_count" -ne 3 ]; then
        echo "ERROR: Expected 3 non-empty rows in '$ped' (found \$row_count)." >&2
        exit 1
    fi

    proband_line=\$(awk '\$3 != "0" && \$4 != "0" { print; exit }' "$ped")
    if [ -z "\$proband_line" ]; then
        echo "ERROR: No row with both parental IDs found in '$ped'." >&2
        exit 1
    fi

    read -r _ SAMPLE_P SAMPLE_F SAMPLE_M _ <<< "\$proband_line"

    cat <<EOF > sample_ids.env
SAMPLE_P=\${SAMPLE_P}
SAMPLE_F=\${SAMPLE_F}
SAMPLE_M=\${SAMPLE_M}
EOF
    """
}
