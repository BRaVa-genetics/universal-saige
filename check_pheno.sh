function is_valid_r_var() {
    pheno=$1

    # Array of characters that are not allowed in R lm variable names (probably overly strict)
    invalid_chars=("+", "-", "*", "/", "^", ":", "~", "(", ")", "[", "]", "$", "@", "!", "%", "#", "&", "=", "?", "|", ";", "<", ">", ",", ".", " ")

    for char in ${invalid_chars[@]}
    do
        if [[ $pheno == *"$char"* ]]
        then
            # If any invalid character is found, print a message and return 1
            echo "The string '$pheno' is not valid for an R variable as it contains the character '$char'."
            return 1
        fi
    done

    # If no invalid characters were found, print a success message and return 0
    echo "The string '$pheno' is valid for an R variable."
    return 0
}
