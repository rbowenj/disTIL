FROM rachelbj/hla-report:1.0.0

RUN Rscript -e 'install.packages(c("DT", "ggplot2"), repos = "http://cran.rstudio.com/");'

ADD ./report_generator.Rmd report_generator.Rmd
#ADD ./report_generator_pdf.Rmd report_generator_pdf.Rmd
ADD ./report_styles.css report_styles.css 

CMD /bin/bash