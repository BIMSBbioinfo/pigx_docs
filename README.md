# Documentation for all PiGx

## Contribute to the development

You can contribute to the development of this guide using GitHub
features such as pull-requests and issue creation.

## Build the book

This book is compiled with pandoc.

Run `guix environment -m manifest.scm` to enter a suitable environment
to hack on this book.  Then run `npm ci` to install the lunr
JavaScript library for the search index.

## How to update the book

Edit the .md files using markdown syntax and run **bash
build.sh**. This will regenerate the book and store the files in the
`book` directory.

Updating the github pages should happen automatically with the `publish` action.
