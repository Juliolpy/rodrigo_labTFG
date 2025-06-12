sphinx-apidoc -o docs src/;
cd docs;

echo "Building html documentation...";
rm -rf _build/html &&  make html;
echo -e "Documentation built successfully!\n\n";
