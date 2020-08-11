mkdir -p tmp

# Build navigation
pandoc -f markdown -t html SUMMARY.md | \
    sed -e 's|\.md|.html|' > book/SUMMARY.html

echo '[' > tmp/documents.json

for file in $(ls -1 pages); do
    echo $file
    title=$(egrep --max-count=1 '^# ' pages/$file | tail -c +3)

    # Build JSON documents for indexing
    tr \" \' < pages/$file | tr -d "[:punct:]" | \
        pandoc -f gfm \
               -t plain \
               -V url="$(basename $file .md).html" \
               -V pagetitle="$title" \
               --template template.json - | tr "\n" " " >> tmp/documents.json

    # Build web pages
    # TODO: level should be taken from SUMMARY
    pandoc -f markdown_github+smart \
           -t html \
           -V level="1.1" \
           -V navigation="$(cat book/SUMMARY.html)" \
           -V pagetitle="$title" \
           -o "book/$(basename $file .md).html" \
           --template template.html pages/$file
done
echo '{}]' >> tmp/documents.json

cp book/README.html book/index.html
node bake-index.js < tmp/documents.json > book/search_index.json
rm tmp/documents.json 
