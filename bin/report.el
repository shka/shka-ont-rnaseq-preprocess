(require 'package)
(package-initialize)
(org-babel-do-load-languages 'org-babel-load-languages '((R . t) (shell . t)))
(setq org-confirm-babel-evaluate nil)
(find-file "out/report.org")
(org-html-export-to-html nil nil nil nil nil)
