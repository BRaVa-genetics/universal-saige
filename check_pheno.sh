# A column name SAIGE can use. It pastes the phenotype and covariate names into
# an R formula without backticks (saige-slim LEDGER #68), so a character with a
# meaning there breaks the fit or silently changes the model ('+' adds a term,
# ':' makes an interaction, '-' removes one). '.' and '_' are fine in R names;
# a name must start with a letter or '.'.
is_valid_r_var () {   # $1 = a column name; prints why and returns 1 if an R formula would misread it
  local name=$1 char
  local invalid_chars=('+' '-' '*' '/' '^' ':' '~' '(' ')' '[' ']' '{' '}' '$' '@' '!' '%' '#' '&' '=' '?' '|' ';' '<' '>' ',' ' ' '"' "'" '`')
  for char in "${invalid_chars[@]}"; do
    if [[ ${name} == *"${char}"* ]]; then
      echo "The column name '${name}' contains '${char}', which an R formula would misread; rename the column."
      return 1
    fi
  done
  if [[ ! ${name} =~ ^[A-Za-z.] ]]; then
    echo "The column name '${name}' must start with a letter or '.' to be an R name; rename the column."
    return 1
  fi
  return 0
}
