/**
 * plugin.js
 *
 * Released under LGPL License.
 * Copyright (c) 1999-2015 Ephox Corp. All rights reserved
 *
 * License: http://www.tinymce.com/license
 * Contributing: http://www.tinymce.com/contributing
 */

tinymce.PluginManager.add("hr",function(n){n.addCommand("InsertHorizontalRule",function(){n.execCommand("mceInsertContent",!1,"<hr />")}),n.addButton("hr",{icon:"hr",tooltip:"Horizontal line",cmd:"InsertHorizontalRule"}),n.addMenuItem("hr",{icon:"hr",text:"Horizontal line",cmd:"InsertHorizontalRule",context:"insert"})});