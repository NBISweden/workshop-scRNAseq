function RawInline (raw)
  if raw.format:match 'html' and raw.text == '<?quarto.version?>'then
    return tostring(quarto.version)
  end
  if raw.format:match 'html' and raw.text == '<?lua.date?>'then
    return tostring(os.date("%d-%m-%Y"))
  end
  if raw.format:match 'html' and raw.text == '<?lua.year?>'then
    return tostring(os.date("%Y"))
  end
  if raw.format:match 'html' and raw.text == '<?lua.time?>'then
    return tostring(os.date("%H:%M:%S"))
  end
end