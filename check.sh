#!/bin/bash

# Set the shared library file
LIBFILE="Rfast.so"

# Ensure the file exists
if [[ ! -f "$LIBFILE" ]]; then
    echo "Error: File $LIBFILE not found!"
    exit 1
fi

echo "Checking if $LIBFILE references R_lsInternal..."
nm -gU "$LIBFILE" | grep R_lsInternal

if [[ $? -eq 0 ]]; then
    echo "✅ Found R_lsInternal in $LIBFILE!"
    exit 0
else
    echo "❌ R_lsInternal not found in $LIBFILE. Checking dependencies..."
fi

# Get the linked libraries
LIBS=$(ldd "$LIBFILE" | awk '{print $3}' | grep -v "(")

if [[ -z "$LIBS" ]]; then
    echo "No dependencies found or ldd failed!"
    exit 1
fi

# Check each linked library for R_lsInternal
for lib in $LIBS; do
    if [[ -f "$lib" ]]; then
        echo "🔍 Searching in $lib..."
        nm -gU "$lib" 2>/dev/null | grep R_lsInternal
        if [[ $? -eq 0 ]]; then
            echo "✅ Found R_lsInternal in $lib!"
            exit 0
        fi
    fi
done

echo "❌ R_lsInternal not found in any linked libraries."
exit 1
