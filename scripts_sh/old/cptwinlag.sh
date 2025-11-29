"""
mkdir -p lagFolder

for f in sigma_*/twin*; do
  [[ -e "$f" ]] || continue                 # salta se non ci sono match
  d="$(dirname "$f")"                       # es. sigma_0.947
  k="${d##*/}"; k="${k#sigma_}"            # -> 0.947 (estrae k)

  base="$(basename "$f")"                   # es. twin.png o twin
  ext="${base##*.}"; name="${base%.*}"

  if [[ "$ext" == "$base" ]]; then          # nessuna estensione
    cp "$f" "lagFolder/${name}_k${k}"
  else
    cp "$f" "lagFolder/${name}_k${k}.${ext}"
  fi
done
"""