# Embedded land geometry

`world_110m.txt` is a compact conversion of Natural Earth's
`ne_110m_land.geojson`, retrieved from the
[`natural-earth-vector` repository](https://github.com/nvkelso/natural-earth-vector/blob/master/geojson/ne_110m_land.geojson)
on 2026-08-08.

Natural Earth [places its vector data in the public
domain](https://www.naturalearthdata.com/about/terms-of-use/). The conversion
keeps polygon and hole order, rounds longitude and latitude to six decimal
places, and drops feature properties that are not used by the land mask. `O`
marks an outer ring and `H` a hole. The renderer embeds this file at compile
time and performs no network request.
