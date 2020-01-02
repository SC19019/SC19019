library(Rcpp) # Attach R package "Rcpp"
# Define function "fc"
cppFunction('int fc(int x, int y, int z) {
int sum = x + y + z;
return sum;
}')