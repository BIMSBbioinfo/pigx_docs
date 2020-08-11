mkdir -p tmp

# Build navigation
pandoc -f markdown -t html -o "book/SUMMARY.html" SUMMARY.md
sed -i book/SUMMARY.html -e 's|\.md|.html|'

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
    pandoc -f markdown_github+smart \
           -t html \
           -V navigation="$(cat book/SUMMARY.html)" \
           -V pagetitle="$title" \
           -o "book/$(basename $file .md).html" \
           --template template.html pages/$file
done
echo '{}]' >> tmp/documents.json

cp book/README.html book/index.html
node bake-index.js < tmp/documents.json > book/search_index.json
rm tmp/documents.json 
