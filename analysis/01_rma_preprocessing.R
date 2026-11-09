# ==========================================
# GSE2034: RMA Preprocessing
# Platform: Affymetrix HG-U133A (GPL96)
# ==========================================


# Load required packages
library(affy)
library(hgu133acdf)

# Define project directories
project_dir <- normalizePath(
  getwd(),
  winslash = "/"
)

cel_dir <- file.path(
  project_dir,
  "data",
  "raw_cel_zipped",
  "cel",
  "extracted"
)

output_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

# Create output directory if necessary
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Find CEL files
cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Check sample count
cat("CEL files found:", length(cel_files), "\n")
stopifnot(length(cel_files) == 286)

# Read raw Affymetrix CEL files
cat("Reading CEL files...\n")

raw_data <- ReadAffy(
  filenames = cel_files,
  cdfname = "HG-U133A"
)

# Verify imported samples
cat("Samples loaded:", length(sampleNames(raw_data)), "\n")
cat("Probes per array:", nrow(exprs(raw_data)), "\n")

# Run RMA preprocessing
cat("Running RMA preprocessing...\n")

normalized_data <- rma(raw_data)

# Extract normalized expression matrix
expression_matrix <- exprs(normalized_data)

# Check output dimensions
cat(
  "RMA matrix dimensions:",
  nrow(expression_matrix), "probe sets x",
  ncol(expression_matrix), "samples\n"
)

# Convert expression matrix to data frame
expression_matrix <- as.data.frame(
  expression_matrix,
  check.names = FALSE
)

# Add probe set identifiers as first column
expression_matrix <- data.frame(
  ID_REF = rownames(expression_matrix),
  expression_matrix,
  check.names = FALSE
)

# Remove row names
rownames(expression_matrix) <- NULL

# Clean sample identifiers
colnames(expression_matrix) <- sub(
  "\\.CEL$",
  "",
  colnames(expression_matrix),
  ignore.case = TRUE
)

# Validate column names
stopifnot(colnames(expression_matrix)[1] == "ID_REF")
stopifnot(ncol(expression_matrix) == 287)
stopifnot(!anyDuplicated(colnames(expression_matrix)))

# Separate numeric expression values for validation
expression_values <- as.matrix(
  expression_matrix[, -1, drop = FALSE]
)

# Validate RMA expression matrix
cat("Missing values:", sum(is.na(expression_values)), "\n")
cat("Non-finite values:", sum(!is.finite(expression_values)), "\n")
cat("Minimum expression:", min(expression_values), "\n")
cat("Maximum expression:", max(expression_values), "\n")
cat("Mean expression:", mean(expression_values), "\n")

# Verify probe identifiers and sample IDs
cat("\nFirst 5 probe IDs:\n")
print(head(expression_matrix$ID_REF, 5))

cat("\nFirst 5 sample IDs:\n")
print(colnames(expression_matrix)[2:6])

# Export compressed RMA expression matrix
output_file <- file.path(
  output_dir,
  "GSE2034_RMA_expression.csv.gz"
)

compressed_file <- gzfile(
  output_file,
  open = "wt"
)

write.csv(
  expression_matrix,
  file = compressed_file,
  quote = FALSE,
  row.names = FALSE
)

close(compressed_file)

cat("\nCompressed RMA file exported successfully:\n")
cat(output_file, "\n")

# Report file size
cat(
  "Compressed file size (MB):",
  round(file.info(output_file)$size / (1024^2), 2),
  "\n"
)
