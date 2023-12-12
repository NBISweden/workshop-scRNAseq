local function ensureHtmlDeps()
  quarto.doc.add_html_dependency({
  name = "reveal-header",
  version = "1.0.0",
  scripts = {
    { path = "resources/js/add_header.js", attribs = {defer = "true"}}
  },
  stylesheets = {"resources/css/add_header.css"}
})
end

if quarto.doc.is_format('revealjs') then
  -- Ensuring the dependencies got loaded before proceeding
  ensureHtmlDeps()
  function Pandoc(doc)
    local blocks = doc.blocks
    local str = pandoc.utils.stringify
    local meta = doc.meta

    -- make divs structure for holding text and logo.
    local header_logo_left = meta['header-logo-left'] and str(meta['header-logo-left']) or ""
    local header_img_left = pandoc.Div(pandoc.Image("", header_logo_left, ""), {class = "header-logo-left"})
    local header_logo_right = meta['header-logo-right'] and str(meta['header-logo-right']) or ""
    local header_img_right = pandoc.Div(pandoc.Image("", header_logo_right, ""), {class = "header-logo-right"})
    local div = pandoc.Div(
      {
        header_img_left,
        header_img_right
      }, 
      {class = 'reveal-header'})
    table.insert(blocks, div)
    return doc
  end
end
