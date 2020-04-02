(require 'package)
(package-initialize)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(setq url-http-attempt-keepalives nil)
(package-refresh-contents)
(package-install 'ess)
(package-install 'htmlize)
(package-install 'yaml-mode)

