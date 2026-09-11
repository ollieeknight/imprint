BEGIN {
  FS = "\t"
  while ((getline line < allowlist) > 0) {
    sub(/\r$/, "", line)
    valid[line] = 1
  }
}
{
  total++
  split($1, q, "_")
  split(q[length(q)], pair, "+")
  db = ""
  for (i = 12; i <= NF; i++) {
    if ($i ~ /^DB:Z:/) {
      db = substr($i, 6)
      break
    }
  }
  split(db, tag, "-")
  if (length(pair) != 2 || length(tag) != 2 || pair[1] != tag[1] || pair[2] != tag[2] || !valid[pair[1]] || !valid[pair[2]]) bad++
}
END {
  print "metric\tvalue"
  print "records\t" total
  print "invalid_or_mismatched_records\t" (bad + 0)
  if (total == 0 || bad > 0) exit 42
}
