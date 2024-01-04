-- quarto revealjs logo extension --

local function ensureHtmlDeps()
  quarto.doc.add_html_dependency({
  name = "reveal-logo",
  version = "1.0.0",
  scripts = {
    { path = "resources/js/reveal-logo.js", attribs = {defer = "true"}}
  },
  stylesheets = {"resources/css/reveal-logo.css"}
})
end

if quarto.doc.is_format('revealjs') then
  -- Ensuring the dependencies got loaded before proceeding
  ensureHtmlDeps()
  function Pandoc(doc)
    local blocks = doc.blocks
    local str = pandoc.utils.stringify
    local meta = doc.meta

    -- read and stringify input if it exists
    local function makeString(path)
      if path == nil or path == '' then
        return nil
      else
        return pandoc.utils.stringify(path)
      end
    end

    -- make image
    local function makeImage(path, height)
      if makeString(height) ~= nil then
        return pandoc.Image("", path, "", pandoc.Attr("", {}, {{"style", "height:" .. height .. ";max-width:none;"}}))
      else
        return pandoc.Image("", path, "")
      end
    end

    -- convert image to link
    local function makeLink(img, url)
      if url ~= nil and url ~= '' then
        return pandoc.Link(img, url)
      end
    end

    -- insert into div
    local function makeDiv(x, class)
      if x == nil or x == '' then
        return x
      else
        return pandoc.Div(x, { class = class })
      end
    end

    -- build a div with logo
    -- @param path Path to an image (string)
    -- @param url A link / url (string)
    -- @param height Height of the image in css units (string)
    -- @param class Class of the div with logo (string)
    --
    local function makeLogo(path, url, height, class)
      if makeString(path) ~= nil and makeString(url) ~= nil then
        return makeDiv(makeLink(makeImage(makeString(path), makeString(height)), makeString(url)), class)
      elseif makeString(path) ~= nil and makeString(url) == nil then
        return makeDiv(makeImage(makeString(path), makeString(height)), class)
      else
        return makeDiv(" ", class)
      end
    end

    -- make divs structure for holding text and logo.
    local header_img_left = makeLogo(meta['header-logo-left'], meta['header-logo-left-url'],
      meta['header-logo-left-height'], "header-logo-left")
    local header_img_right = makeLogo(meta['header-logo-right'], meta['header-logo-right-url'],
      meta['header-logo-right-height'], "header-logo-right")

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
