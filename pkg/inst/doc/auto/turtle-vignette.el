(TeX-add-style-hook "turtle-vignette"
 (lambda ()
    (LaTeX-add-bibliographies
     "turtle")
    (LaTeX-add-labels
     "fig:data1"
     "fig:condensed"
     "fig:cml1"
     "sec:quickstart")
    (TeX-add-symbols
     '("fixme" 1)
     "R"
     "Splus")
    (TeX-run-style-hooks
     "alltt"
     "url"
     "babel"
     "american"
     "palatino"
     "latex2e"
     "art11"
     "article"
     "11pt")))

