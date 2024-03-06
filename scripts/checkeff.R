chekEff <- function(data, ..., x, y, fill = NULL, size = 1, scales = "fixed") {
  
  # check for non-continuous x, group if necessary
  if(is.character(data %>% dplyr::pull({{x}})) |
     is.factor(data %>% dplyr::pull({{x}})) |
     is.logical(data %>% dplyr::pull({{x}}))) {
    p <- ggplot2::ggplot(data) +
      ggplot2::aes(x = {{x}}, y = {{y}}) +
      ggplot2::facet_wrap(ggplot2::vars(...), scales = scales) +
      theme_bw() +
      theme(strip.background = ggplot2::element_rect(
        color="black", fill="white", size=0.5, linetype="solid"),
        panel.grid = element_blank(),
        strip.text = ggplot2::element_text(face = "italic", size = 5))  # Make strip text italic
    
  } else {
    p <- ggplot2::ggplot(data) +
      ggplot2::aes(x = {{x}}, y = {{y}}, group = {{x}}) +
      ggplot2::facet_wrap(ggplot2::vars(...), scales = scales) +
      theme_bw() +
      theme(strip.background = ggplot2::element_rect(
        color="black", fill="white", size=0.5, linetype="solid"),
        panel.grid = element_blank(),
        strip.text = ggplot2::element_text(face = "italic", size = 5))  # Make strip text italic
  }
  
 

  # check on optional fill parameter
  if(as.character(rlang::enquo(fill))[2] != "NULL"){
    p.out <- p + ggplot2::geom_jitter(aes(fill = as.character({{fill}})),
                                      position = ggplot2::position_jitter(),
                                      size = size,
                                      shape = 21) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, fill = "white") +
      ggplot2::labs(fill = rlang::enquo(fill)) +
      ggplot2::scale_x_continuous(breaks = c(0, 30), labels = c(0, 30)) +  # Set x-axis breaks and labels
      ggplot2::xlab("Concentration") +  # Set x-axis label
      ggplot2::ylab("Median animal length") +  # Set y-axis label
      ggplot2::theme(legend.position = "none") +
      theme(
        strip.text = ggplot2::element_text(face = "italic"))  # Remove legend
  } else {
    p.out <- p + ggplot2::geom_jitter(width = 0.25,
                                      size = size,
                                      fill = "grey",
                                      shape = 21) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, fill = "white" ) +
      ggplot2::scale_x_continuous(breaks = c(0, 30), labels = c(0, 30)) +  # Set x-axis breaks and labels
      ggplot2::xlab("Concentration") +  # Set x-axis label
      ggplot2::ylab("Median animal length") +  # Set y-axis label
      ggplot2::theme(legend.position = "none")+
      theme(
        strip.text = ggplot2::element_text(face = "italic"))  # Remove legend
  }

  # return it
  return(p.out)
}
