# blosum
BLOSUM generator with the addition of any generic characters for special amino
acids

Example usage:

```blosum --identity 62 --remove-characters "-." --all-remove-characters --field 1 --csv "ARNDCQEGHILKMFPSTWYV" --input blocks.fasta```

This command takes in the blocks database formatted in a fasta file where the
first field of the header contains the block that sequence is in. All sequences
were removed here if they contained at least one of the characters "." or "-".
The output BLOSUM takes into account all characters, even characters that might
represent something special (like ":"). This command outputs a csv file as a 2D
formatted BLOSUM, where the ordering is "ARNDCQEGHILKMFPSTWYV".
