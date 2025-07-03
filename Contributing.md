# Contributing to SpatialMetabolics

We welcome contributions to the SpatialMetabolics package! This document provides guidelines for contributing.

## How to Contribute

### Reporting Issues

- Use the [GitHub issue tracker](https://github.com/yourusername/SpatialMetabolics/issues)
- Check if the issue already exists
- Include a minimal reproducible example
- Provide session info (`sessionInfo()`)
- Describe expected vs actual behavior

### Suggesting Features

- Open an issue with the "enhancement" label
- Describe the feature and its use case
- Provide examples if possible

### Contributing Code

1. **Fork the repository**
   ```bash
   git clone https://github.com/yourusername/SpatialMetabolics.git
   cd SpatialMetabolics
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make changes**
   - Follow the coding style guide
   - Add tests for new features
   - Update documentation
   - Add yourself to `inst/AUTHORS` if substantial contribution

4. **Test your changes**
   ```bash
   make test
   make check
   ```

5. **Submit a pull request**
   - Describe your changes
   - Reference any related issues
   - Ensure all checks pass

## Coding Standards

### R Code Style

We follow the [Bioconductor style guide](https://bioconductor.org/developers/how-to/coding-style/):

- Use 4 spaces for indentation (no tabs)
- Maximum line length: 80 characters
- Use `<-` for assignment, not `=`
- Function names: `camelCase` for exported, `.lowerCamelCase` for internal
- Variable names: `snake_case`

Example:
```r
#' Calculate metabolic scores
#' 
#' @param object SpatialMetabolic object
#' @param pathways List of pathways
#' @return Modified object
calculateMetabolicScores <- function(object, pathways) {
    # Check inputs
    if (!is(object, "SpatialMetabolic")) {
        stop("object must be a SpatialMetabolic object")
    }
    
    # Process pathways
    pathway_scores <- .calculateScores(object, pathways)
    
    # Store results
    metabolicScores(object) <- pathway_scores
    
    return(object)
}
```

### Documentation

- All exported functions must have roxygen2 documentation
- Include `@examples` for user-facing functions
- Use `@keywords internal` for non-exported functions
- Update vignettes for major features

### Testing

- Write unit tests for all new functions
- Aim for >80% code coverage
- Test edge cases and error conditions
- Use `testthat` framework

Example test:
```r
test_that("calculateMetabolicScores works correctly", {
    # Setup
    spm <- createTestData()
    
    # Test basic functionality
    spm_scores <- calculateMetabolicScores(spm, pathways = list(test = c("Gene1", "Gene2")))
    expect_s4_class(spm_scores, "SpatialMetabolic")
    expect_true(nrow(metabolicScores(spm_scores)) > 0)
    
    # Test error handling
    expect_error(calculateMetabolicScores("not an object"), "must be a SpatialMetabolic")
})
```

## Development Workflow

### Setting Up Development Environment

1. Install development dependencies:
   ```r
   install.packages(c("devtools", "testthat", "roxygen2", "BiocCheck"))
   BiocManager::install(c("BiocStyle", "BiocGenerics"))
   ```

2. Load package for development:
   ```r
   devtools::load_all()
   ```

### Common Tasks

- **Add a new function**: Create in appropriate R/ file, add tests, document
- **Fix a bug**: Add failing test first, then fix
- **Update documentation**: Edit .Rd files via roxygen2, update vignettes
- **Add a dataset**: Place in inst/extdata/ with documentation

### Before Submitting

1. Run full check:
   ```bash
   make check
   ```

2. Run BiocCheck:
   ```bash
   make bioccheck
   ```

3. Update NEWS.md with your changes

4. Ensure all tests pass:
   ```bash
   make test
   ```

## Getting Help

- Read existing code for style examples
- Check the [Bioconductor developer guide](https://bioconductor.org/developers/)
- Ask questions in issues or discussions
- Email: your.email@example.com

## Code of Conduct

### Our Pledge

We pledge to make participation in our project a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4.

## Recognition

Contributors will be acknowledged in:
- The package DESCRIPTION file
- The inst/AUTHORS file
- Release announcements

Thank you for contributing to SpatialMetabolics!
