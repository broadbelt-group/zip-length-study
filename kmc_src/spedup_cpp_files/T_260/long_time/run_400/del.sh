for f in chain_weights_*_mins.csv; do
    n=$(echo "$f" | sed -E 's/chain_weights_([0-9]+)_mins\.csv/\1/')
    if (( n % 100 != 0 )); then
        echo "Deleting $f"
        rm "$f"
    fi
done

