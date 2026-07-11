# Design Descriptions for DoEIY

get_design_description <- function(design_type) {
  descs <- list(
    "Plackett-Burman" = '<p style="text-align: justify;"><strong>Plackett-Burman Designs</strong> are efficient screening tools for estimating main effects in an experiment. These designs allow for independent estimation of main effects, though they may be confounded with higher-order interactions. This software offers 2-level Plackett-Burman designs suitable for experiments involving 4 to 23 factors.</p>',
    "Full Factorial" = '<p style="text-align: justify;"><strong>Full Factorial Designs</strong> investigate every possible combination of factor levels, using two or more levels for each factor. This comprehensive approach enables the independent estimation of all main effects and higher-order interactions, providing a complete understanding of factor relationships.</p>',
    "Fractional Factorial" = "<p style='text-align: justify;'><strong>Fractional Factorial Designs</strong> are efficient alternatives to full factorial designs, utilizing only a fraction of the total possible runs. They allow for the estimation of main effects and, depending on the design's resolution, certain interactions. This software provides 2-level fractional factorial designs with the following resolutions: </p>
          <ul>
          <li style='text-align: justify;'><strong>Resolution III</strong> - main effects are not confounded with each other but may be confounded with two-factor interactions.</li>
          <li style='text-align: justify;'><strong>Resolution IV</strong> - main effects are not confounded with other main effects or two-factor interactions; two-factor interactions may still be confounded with each other.</li>
          <li style='text-align: justify;'><strong>Resolution V</strong> - two-factor interactions are not confounded with main effects or other two-factor interactions, but may be confounded with three-factor interactions.</li>
          </ul>",
    "Box-Behnken" = "<p><strong>Box-Behnken Designs</strong> are used to estimate quadratic models by employing three levels for each factor. This design enables the independent estimation of main effects, two-factor interactions, and quadratic effects. This software offers Box-Behnken designs for experiments with 3 to 7 factors.</p>",
    "Central Composite" = "<p style='text-align: justify;'><strong>Central Composite Designs (CCD)</strong> are response surface designs used to estimate quadratic models. They consist of a factorial or fractional factorial design, additional 'star' (axial) points, and center points. This structure allows for the estimation of curvature and the fitting of second-order response surfaces. This software offers Circumscribed, Inscribed, and Face Centered central composite designs for experiments with 2 or more factors.</p>",
    "Latin Hypercube Sampling" = "<p style='text-align: justify;'><strong>Latin Hypercube Sampling (LHS)</strong> is a space-filling design used to explore multi-dimensional parameter spaces. Unlike traditional factorial designs, LHS ensures that each factor is uniformly sampled across its range by dividing the range of each factor into equal intervals and selecting one point from each interval. This software provides improved LHS designs, optimizing the space-filling properties, suitable for computer experiments and simulations with any number of factors.</p>",
    "D-Optimal" = "<p style='text-align: justify;'><strong>D-Optimal Designs</strong> are computer-generated designs optimized based on a specified model. They choose a subset of runs from a candidate set of all possible combinations to minimize the variance of the estimated model coefficients (specifically maximizing the determinant of the information matrix). This approach is highly flexible, allowing for custom model specifications, constraints on run sizes, and the inclusion of continuous, discrete, and categorical factors with unequal levels.</p>"
  )
  
  if (design_type %in% names(descs)) {
    return(descs[[design_type]])
  }
  return("")
}
