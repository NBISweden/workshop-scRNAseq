---
layout: default
title: Single cell RNA sequencing analysis course
---

# Single cell RNA sequencing analysis course
Scilifelab Solna, Rooms Air & Fire,  2018-05-21 - 2018-05-23

{% highlight bash %}
# check alignment of text
asbd      ösjgödlkj
ögkjlök   jgkldj
kjk       llll
ii        asg
{% endhighlight %}

``` r
# test R code block
nPlot <- 5 #number of genes top plot per pc.
par(mfrow=c(5,2),mar=c(2,2,1,1),oma=c(1,5,1,1))
for (i in 1:5) {
    top<-order(PC$rotation[,i],decreasing=T)[1:nPlot]
    bottom<-order(PC$rotation[,i],decreasing=F)[1:nPlot]
    barplot(contr[top,i],main=sprintf("genes on pos axis PC%d",i),
            ylab="% contr",las=2,horiz=T)
    barplot(contr[bottom,i],main=sprintf("genes on neg axis PC%d",i),
            ylab="% contr",las=2,horiz=T)
}
```

{% highlight bash %}
# Johan code block with 
$ bundle exec jekyll serve
mmmmm mmmmm
lllll lllll
{% endhighlight %}


```js
// Javascript code with syntax highlighting.
var fun = function lang(l) {
  dateformat.i18n = require('./lang/' + l)
  return true;
}
```

```ruby
# Ruby code with syntax highlighting
GitHubPages::Dependencies.gems.each do |gem, version|
  s.add_dependency(gem, "= #{version}")
end
```



##### Schedule

Course schedule can be found here: [Schedule](schedule)

##### Exercises

All exercises for the afternoon sessions can be found at [Exercises](exercises). 

For working on Uppmax: to use the allocations we have for the course, please look [here](login.md).

##### Precourse

Please read carefully the [Precourse material](precourse) before the course start. 

##### Address and travel suggestions

[Travel Info](travel)


##### Course leaders

* [Åsa Björklund](http://nbis.se/about/staff/asa-bjorklund/)
* [Stefania Giacomello](http://nbis.se/about/staff/stefania-giacomello/)

