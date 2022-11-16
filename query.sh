# Search API
curl -i -H "Accept: application/json" -H "Content-Type: application/json" -X GET https://api.gbif.org/v1/species/7707728/vernacularNames
curl -i -H "Accept: application/json" -H "Content-Type: application/json" -X GET https://api.gbif.org/v1/species/match?name=Trochilidae

# Download API
curl --include --user username:password --header "Content-Type: application/json" --data @query.json https://api.gbif.org/v1/occurrence/download/request


