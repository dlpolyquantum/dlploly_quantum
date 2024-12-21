#!/bin/bash

# Define lists
list1=("C" "N" "O" "F" "S" "Li")
list2=("OW" "HW")

nrows=998
step=0.01
# Generate 4-column data (columns are index, scaled value, 0.00, 0.00)
# Save this to a temporary data file to be reused
temp_data="temp_data.txt"
for i in $(seq 1 $nrows); do
    value=$(echo "scale=2; $i * $step" | bc)
    printf "%-4d %12.6f %12.6f %12.6f\n" "$i" "$value" 0.00 0.00 >> "$temp_data"
done

# Counter for naming files
counter=1

# Loop through combinations of list1 with list2 and list2 with itself
for item1 in "${list1[@]}"; do
    for item2 in "${list2[@]}"; do
        # Create the filename
        filename="Pair_${counter}_pot.table"
        
        # Write the pair header and data to the file
        echo "$item1 $item2" > "$filename"
        cat "$temp_data" >> "$filename"
        
        # Increment counter
        ((counter++))
    done
done

# Loop through combinations within list2 itself
for i in "${!list2[@]}"; do
    for j in $(seq "$i" $((${#list2[@]} - 1))); do
        # Create the filename
        filename="Pair_${counter}_pot.table"
        
        # Write the pair header and data to the file
        echo "${list2[i]} ${list2[j]}" > "$filename"
        cat "$temp_data" >> "$filename"
        
        # Increment counter
        ((counter++))
    done
done

# Clean up temporary data file
rm "$temp_data"

echo "Files generated: Pair_1_pot.table to Pair_$((counter - 1))_pot.table"


# Loop to generate Potential.dat with zero values
rm -f "Potential.dat"

for i in $(seq 0 $nrows); do
    index=$(awk "BEGIN {printf \"%.2f\", $i * $step }")
    # Create a line with the index followed by number of zeros
    echo "$index $(printf ' 0.000 %.0s' $(seq 2 $counter))" >> "Potential.dat"
done

