pandoc -f markdown -t html -o "book/SUMMARY.html" SUMMARY.md
sed -i book/SUMMARY.html -e 's|\.md|.html|'

for file in $(ls -1 pages); do
    echo $file
    pandoc -f markdown_github+smart \
           -t html \
           -V navigation="$(cat book/SUMMARY.html)" \
           -o "book/$(basename $file .md).html" \
           --template template.html pages/$file
done

cp book/README.html book/index.html

