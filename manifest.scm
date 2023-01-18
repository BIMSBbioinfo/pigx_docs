;; Use guix environment -m manifest.scm
(specifications->manifest
 (list "pandoc"
       "node"
       "bash"
       "coreutils"
       "sed"
       "rsync"))
