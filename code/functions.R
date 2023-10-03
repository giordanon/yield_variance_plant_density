
s1 <- rel(1.5) 
theme_custom <-  function(size = rel(1.5)){


theme <- theme(axis.title.y = element_text(size = s1, angle = 90), 
              axis.title.x = element_text(size = s1), 
              axis.ticks = element_line(size = 0.2),
              axis.ticks.length = unit(0.15, "cm"),
              panel.background = element_rect(fill = "white",
                                              colour = "black",
                                              size = 0.5, linetype = "solid"),
              strip.text.x = element_text(size = s1),
              strip.placement = "outside",
              axis.text.y = element_text(hjust = 1, size = s1), 
              axis.text.x = element_text(hjust = 0.5, angle = 0,  size = s1),
              strip.background =element_rect(fill="white"),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.title=element_text(size = s1),
              legend.text=element_text(size = s1), 
              legend.background = element_rect(fill = "white", color = "black"), 
              )
return(theme)
}


mod_predict <- function(posterior){
  
  a1 <- posterior[,c("a1")]
  maxVPD <- posterior[,c("maxVPD")]
  minVPD <- posterior[,c("minVPD")]
  
  b1 <- posterior[,"b1"]
  AOPD <- posterior[,"AOPD"]
  
  # Compare yields and variance among treatments at the expected value of AOPD and MRPD
  q500 <- apply(posterior , MARGIN = 2, FUN = mean)
  
  var <- expPlateau(a1, maxVPD, minVPD, rep(q500["minVPD"], nrow(posterior)) )
  
  mean <- quadPlateau(b1, AOPD, rep(q500["AOPD"], nrow(posterior)))
  
  sd <- sqrt(var)
  
  return(cbind("sd" = sd,"mean" = mean))
  
}





