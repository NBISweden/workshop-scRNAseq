-- Import the os and io modules
os = require('os')
io = require('io')

-- The filter function
function Meta(meta)
  -- Execute 'quarto version' command and get the output
  local quarto_version = io.popen('quarto --version'):read('*all')
  if quarto_version then
    -- Add the version to the document metadata
    meta['quarto_version'] = quarto_version
  end

  -- Get the current date and time
  local current_date = os.date("%d-%m-%Y")
  local current_year = os.date("%Y")
  local current_time = os.date("%H:%M:%S")

  -- Add them to the document metadata
  meta['current_date'] = current_date
  meta['current_year'] = current_year
  meta['current_time'] = current_time

  return meta
end
