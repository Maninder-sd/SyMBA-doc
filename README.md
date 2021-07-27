# SyMBA-doc

## things to do if doc is updated
* make changes in google docs (remove margins in doc)
* download the html of the doc (file>download)
* place into doc-html
The html is displayed using iframe on the actual webpage (embedding the html into the website)
* change iframe src (on actual website) if needed
* the font being used in Open Sans, need to import it into the downloaded html
  * open generated html, select all, press ctrl+shift+i to format it 
  *  The first couple line should look like this:

\<html>

\<head>

    <meta content="text/html; charset=UTF-8" http-equiv="content-type">
    
    <style type="text/css">
    
    @import url('https://fonts.googleapis.com/css2?family=Open+Sans:ital,wght@0,300;0,400;0,600;0,700;1,400&display=swap');
    
  * If not, copy-paste the import right after the <style> tag

@import url('https://fonts.googleapis.com/css2?family=Open+Sans:ital,wght@0,300;0,400;0,600;0,700;1,400&display=swap');

* Any hyperlink in the generated doc needs the target attribute
  * ctrl+f to search '<a', this highlights all areas where hyperlinks are 
  * copy-paste target="_blank" into the element. Should look something like this:
  
\<a target="_blank"

class="c14"

...


