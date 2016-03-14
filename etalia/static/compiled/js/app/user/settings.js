define("app/util/utils",["jquery"],function(e){function t(e){var t=null;if(document.cookie&&""!=document.cookie)for(var i=document.cookie.split(";"),n=0;n<i.length;n++){var s=jQuery.trim(i[n]);if(s.substring(0,e.length+1)==e+"="){t=decodeURIComponent(s.substring(e.length+1));break}}return t}function i(e){return/^(GET|HEAD|OPTIONS|TRACE)$/.test(e)}var n=t("csrftoken");return e.ajaxSetup({beforeSend:function(e,t){i(t.type)||this.crossDomain||e.setRequestHeader("X-CSRFToken",n)}}),{toggleClass:function(e,t,i){return"undefined"!=typeof i?(i?e.addClass(t):e.removeClass(t),i):e.hasClass(t)?(e.removeClass(t),!1):(e.addClass(t),!0)},getParameterByName:function(e,t){t||(t=window.location.href),e=e.replace(/[\[\]]/g,"\\$&");var i=new RegExp("[?&]"+e+"(=([^&#]*)|&|#|$)"),n=i.exec(t);return n?n[2]?decodeURIComponent(n[2].replace(/\+/g," ")):"":null},bindLoadingButtons:function(t){t.on("click",".loading-button",function(t){var i=e(t.target).closest(".loading-button");return i.data("busy")?(t.stopPropagation(),void t.stopImmediatePropagation()):(i.data("busy",!0),void i.find(".eai").removeAttr("class").addClass("eai eai-loading"))})},restoreLoadingButton:function(e,t){e.data("busy",!1),e.find(".eai").removeAttr("class").addClass("eai "+t)},popup:function(e,t,i,n){t=t||"etalia-popup",i=i||520,n=n||460;var s=Math.round(window.innerWidth/2-i/2),a=Math.round(window.innerHeight/2-n/2),o="menubar=no,toolbar=no,resizable=yes,scrollbars=yes,width="+i+",height="+n+",top="+a+",left="+s;return window.open(e,t,o)}}}),define("app/ui/layout/flap",["jquery"],function(e){var t=function(t){t=t||{},this.config=e.extend({debug:!1,flap:null,side:null,button:null,backdrop:null,mobileMaxWidth:992},t),this.$flap=e(this.config.flap),this.$button=e(this.config.button),this.$backdrop=e(this.config.backdrop),this.mobile=null,this.clickHandlersEnabled=!1,this.opened=!1;var i=this;this.log=function(e){i.config.debug&&console.log("[Flap "+this.config.side+"] "+e)},this.buttonClickHandler=function(){i.log("Button click handler"),i.open()},this.backdropClickHandler=function(){i.log("Overlay click handler"),i.close()},this.$flap.on("redraw",function(){i.init()})};return t.prototype.init=function(){this.log("Flap init"),this.$flap.css({height:window.innerHeight+"px"});var e=window.innerWidth<this.config.mobileMaxWidth;e!=this.mobile&&(this.close(),e?this.enableClickHandlers():this.disableClickHandlers()),this.mobile=e},t.prototype.enableClickHandlers=function(){this.clickHandlersEnabled||(this.log("enableClickHandlers"),this.$button.on("click",this.buttonClickHandler),this.$backdrop.on("click",this.backdropClickHandler),this.clickHandlersEnabled=!0)},t.prototype.disableClickHandlers=function(){this.clickHandlersEnabled&&(this.log("disableClickHandlers"),this.$button.off("click",this.buttonClickHandler),this.$backdrop.off("click",this.backdropClickHandler),this.clickHandlersEnabled=!1)},t.prototype.open=function(){return this.opened?this:(e("body").addClass("flap-opened").addClass(this.config.side),void(this.opened=!0))},t.prototype.close=function(){return this.opened?(e("body").removeClass("flap-opened").removeClass(this.config.side),void(this.opened=!1)):this},t}),define("app/ui/layout/interact",["jquery"],function(e){var t=function(t){this.config=e.extend({debug:!1},t)};return t.prototype.log=function(){return this.config.debug&&console.log("[Interact] "+arguments[0],Array.prototype.splice.call(arguments,1)),this},t.prototype.init=function(){this.log("init()");var e=this,t=0,i=setInterval(function(){return e.log("check ga",t,typeof ga),"function"==typeof ga?(e.initHandlers(),void clearInterval(i)):(t++,void(t>20&&clearInterval(i)))},500);return this},t.prototype.initHandlers=function(){this.log("initHandlers()");var t=this,i=e("body"),n="control";return i.on("etalia.control.search.change",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:"search",eventLabel:i.expression};t.log("ga.send()",s),ga("send",s)}).on("etalia.control.timespan.change",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:"timespan",eventLabel:i.label,eventValue:i.value};t.log("ga.send()",s),ga("send",s)}).on("etalia.control.cluster.change",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:"cluster",eventLabel:i.label,eventValue:i.value};t.log("ga.send()",s),ga("send",s)}).on("etalia.control.pinned.change",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:i.active?"pin":"unpin"};t.log("ga.send()",s),ga("send",s)}).on("etalia.control.filters.change",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:(i.active?"select_":"deselect_")+i.group,eventLabel:i.label,eventValue:i.value};t.log("ga.send()",s),ga("send",s)}).on("etalia.control.filters.more",function(e,i){var s={hitType:"event",eventCategory:n,eventAction:"more_filter",eventLabel:i.group};t.log("ga.send()",s),ga("send",s)}),i.on("etalia.publication.pin",function(e,i){var n={hitType:"event",eventCategory:i.getExtra("source"),eventAction:i.isPinned()?"pin":"unpin",eventLabel:i.getExtra("title"),eventValue:i.getId()};t.log("ga.send()",n),ga("send",n)}).on("etalia.publication.ban",function(e,i){var n={hitType:"event",eventCategory:i.getExtra("source"),eventAction:"ban",eventLabel:i.getExtra("title"),eventValue:i.getId()};t.log("ga.send()",n),ga("send",n)}).on("etalia.publication.add",function(e,i){var n={hitType:"event",eventCategory:i.getExtra("source"),eventAction:"add",eventLabel:i.getExtra("title"),eventValue:i.getId()};t.log("ga.send()",n),ga("send",n)}).on("etalia.publication.trash",function(e,i){var n={hitType:"event",eventCategory:i.getExtra("source"),eventAction:"trash",eventLabel:i.getExtra("title"),eventValue:i.getId()};t.log("ga.send()",n),ga("send",n)}).on("etalia.publication.restore",function(e,i){var n={hitType:"event",eventCategory:i.getExtra("source"),eventAction:"restore",eventLabel:i.getExtra("title"),eventValue:i.getId()};t.log("ga.send()",n),ga("send",n)}).on("etalia.publication.trash-clear",function(){var e={hitType:"event",eventCategory:"library",eventAction:"trash-clear"};t.log("ga.send()",e),ga("send",e)}).on("etalia.publication.share",function(e,i){var n={hitType:"event",eventCategory:i.source,eventAction:"share_"+i.support,eventLabel:i.title,eventValue:i.id};t.log("ga.send()",n),ga("send",n)}),i.on("etalia.detail.similar_timespan.change",function(e,i){var n={hitType:"event",eventCategory:"detail",eventAction:"similar_timespan",eventLabel:i.title,eventValue:i.value};t.log("ga.send()",n),ga("send",n)}),this},(new t).init()}),define("app/ui/layout/invite",["jquery"],function(e){function t(t){t.preventDefault();var n=i.find(".modal-form"),s=i.find(".modal-result");return e.ajax({type:n.attr("method"),url:n.attr("action"),data:n.serialize(),success:function(e){i.find(".form-errors").empty(),n.hide().find("input, textarea").val(""),s.show()},error:function(t){console.log(t.responseText);var n=i.find(".form-errors").empty(),s=JSON.parse(t.responseText);e.each(s,function(e,t){n.prepend('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>'+t+"</div>")})}}),!1}var i;e(function(){i=e("#invite-modal"),e("#open-invite").on("click",function(){e("#profile-dropdown").hide(),i.modal("show")}),i.find("form[data-async]").on("submit",t),i.find(".form-errors").empty(),i.on("hidden.bs.modal",function(){i.find(".modal-form").show(),i.find(".modal-result").hide()})})}),define("app/ui/layout/close-alerts",["jquery"],function(e){e(function(){e("a.close[close-href]").click(function(t){t.preventDefault(),e.post(e(this).attr("close-href"),"",function(){})})})}),define("app/ui/layout",["jquery","app/util/utils","app/ui/layout/flap","app/ui/layout/interact","app/ui/layout/invite","app/ui/layout/close-alerts"],function(e,t,i){var n=function(t){this.config=e.extend({debug:!1,leftFlap:"#nav-flap",leftFlapButton:"#toggle-nav",rightFlap:"#filter-flap",rightFlapButton:"#toggle-filter",backdrop:"#backdrop"},t),this.leftFlap=null,this.rightFlap=null,this.busy=!1;var i=this;this.log=function(e){i.config.debug&&console.log("[Layout] "+e)},this.resizeTimeout=null,this.resizeHandler=function(){i.log("Resize handler"),i.leftFlap&&i.leftFlap.init(),i.rightFlap&&i.rightFlap.init()}};n.prototype.init=function(){var n=e(this.config.leftFlap),s=e(this.config.leftFlapButton),a=e(this.config.rightFlap),o=e(this.config.rightFlapButton),r=!1;if(this.log("Init"),n.length?(this.log("Left flap found"),this.leftFlap=new i({debug:this.config.debug,flap:this.config.leftFlap,side:"left",button:this.config.leftFlapButton,backdrop:this.config.backdrop,mobileMaxWidth:992}),r=!0):s.hide(),a.length?(this.log("Right flap found"),this.rightFlap=new i({debug:this.config.debug,flap:this.config.rightFlap,side:"right",button:this.config.rightFlapButton,backdrop:this.config.backdrop,mobileMaxWidth:1200}),r=!0):o.hide(),r){var l=this;e(window).on("resize",function(){l.resizeTimeout&&clearTimeout(l.resizeTimeout),l.resizeTimeout=setTimeout(l.resizeHandler,100)}),l.resizeHandler()}var u=e("#toggle-profile"),h=e("#profile-dropdown");h.length&&(u.on("click",function(i){t.toggleClass(e(i.delegateTarget),"active")?h.show():h.hide()}),h.on("click",function(e){e.stopPropagation()}),e(window).on("click",function(t){0==e(t.target).closest("#toggle-profile").length&&(e("#toggle-profile").removeClass("active"),e("#profile-dropdown").hide())}))},n.prototype.setBusy=function(t){e("#busy-content").html(t),this.busy||(this.busy=!0,e("body").addClass("busy"))},n.prototype.setAvailable=function(){this.busy&&(this.busy=!1,e("#busy-content").html(""),e("body").removeClass("busy"))};var s=new n;return s.init(),s}),function(e){"function"==typeof define&&define.amd?define("jquery-ui/core",["jquery"],e):e(jQuery)}(function(e){function t(t,n){var s,a,o,r=t.nodeName.toLowerCase();return"area"===r?(s=t.parentNode,a=s.name,t.href&&a&&"map"===s.nodeName.toLowerCase()?(o=e("img[usemap='#"+a+"']")[0],!!o&&i(o)):!1):(/^(input|select|textarea|button|object)$/.test(r)?!t.disabled:"a"===r?t.href||n:n)&&i(t)}function i(t){return e.expr.filters.visible(t)&&!e(t).parents().addBack().filter(function(){return"hidden"===e.css(this,"visibility")}).length}e.ui=e.ui||{},e.extend(e.ui,{version:"1.11.4",keyCode:{BACKSPACE:8,COMMA:188,DELETE:46,DOWN:40,END:35,ENTER:13,ESCAPE:27,HOME:36,LEFT:37,PAGE_DOWN:34,PAGE_UP:33,PERIOD:190,RIGHT:39,SPACE:32,TAB:9,UP:38}}),e.fn.extend({scrollParent:function(t){var i=this.css("position"),n="absolute"===i,s=t?/(auto|scroll|hidden)/:/(auto|scroll)/,a=this.parents().filter(function(){var t=e(this);return n&&"static"===t.css("position")?!1:s.test(t.css("overflow")+t.css("overflow-y")+t.css("overflow-x"))}).eq(0);return"fixed"!==i&&a.length?a:e(this[0].ownerDocument||document)},uniqueId:function(){var e=0;return function(){return this.each(function(){this.id||(this.id="ui-id-"+ ++e)})}}(),removeUniqueId:function(){return this.each(function(){/^ui-id-\d+$/.test(this.id)&&e(this).removeAttr("id")})}}),e.extend(e.expr[":"],{data:e.expr.createPseudo?e.expr.createPseudo(function(t){return function(i){return!!e.data(i,t)}}):function(t,i,n){return!!e.data(t,n[3])},focusable:function(i){return t(i,!isNaN(e.attr(i,"tabindex")))},tabbable:function(i){var n=e.attr(i,"tabindex"),s=isNaN(n);return(s||n>=0)&&t(i,!s)}}),e("<a>").outerWidth(1).jquery||e.each(["Width","Height"],function(t,i){function n(t,i,n,a){return e.each(s,function(){i-=parseFloat(e.css(t,"padding"+this))||0,n&&(i-=parseFloat(e.css(t,"border"+this+"Width"))||0),a&&(i-=parseFloat(e.css(t,"margin"+this))||0)}),i}var s="Width"===i?["Left","Right"]:["Top","Bottom"],a=i.toLowerCase(),o={innerWidth:e.fn.innerWidth,innerHeight:e.fn.innerHeight,outerWidth:e.fn.outerWidth,outerHeight:e.fn.outerHeight};e.fn["inner"+i]=function(t){return void 0===t?o["inner"+i].call(this):this.each(function(){e(this).css(a,n(this,t)+"px")})},e.fn["outer"+i]=function(t,s){return"number"!=typeof t?o["outer"+i].call(this,t):this.each(function(){e(this).css(a,n(this,t,!0,s)+"px")})}}),e.fn.addBack||(e.fn.addBack=function(e){return this.add(null==e?this.prevObject:this.prevObject.filter(e))}),e("<a>").data("a-b","a").removeData("a-b").data("a-b")&&(e.fn.removeData=function(t){return function(i){return arguments.length?t.call(this,e.camelCase(i)):t.call(this)}}(e.fn.removeData)),e.ui.ie=!!/msie [\w.]+/.exec(navigator.userAgent.toLowerCase()),e.fn.extend({focus:function(t){return function(i,n){return"number"==typeof i?this.each(function(){var t=this;setTimeout(function(){e(t).focus(),n&&n.call(t)},i)}):t.apply(this,arguments)}}(e.fn.focus),disableSelection:function(){var e="onselectstart"in document.createElement("div")?"selectstart":"mousedown";return function(){return this.bind(e+".ui-disableSelection",function(e){e.preventDefault()})}}(),enableSelection:function(){return this.unbind(".ui-disableSelection")},zIndex:function(t){if(void 0!==t)return this.css("zIndex",t);if(this.length)for(var i,n,s=e(this[0]);s.length&&s[0]!==document;){if(i=s.css("position"),("absolute"===i||"relative"===i||"fixed"===i)&&(n=parseInt(s.css("zIndex"),10),!isNaN(n)&&0!==n))return n;s=s.parent()}return 0}}),e.ui.plugin={add:function(t,i,n){var s,a=e.ui[t].prototype;for(s in n)a.plugins[s]=a.plugins[s]||[],a.plugins[s].push([i,n[s]])},call:function(e,t,i,n){var s,a=e.plugins[t];if(a&&(n||e.element[0].parentNode&&11!==e.element[0].parentNode.nodeType))for(s=0;s<a.length;s++)e.options[a[s][0]]&&a[s][1].apply(e.element,i)}}}),function(e){"function"==typeof define&&define.amd?define("jquery-ui/widget",["jquery"],e):e(jQuery)}(function(e){var t=0,i=Array.prototype.slice;return e.cleanData=function(t){return function(i){var n,s,a;for(a=0;null!=(s=i[a]);a++)try{n=e._data(s,"events"),n&&n.remove&&e(s).triggerHandler("remove")}catch(o){}t(i)}}(e.cleanData),e.widget=function(t,i,n){var s,a,o,r,l={},u=t.split(".")[0];return t=t.split(".")[1],s=u+"-"+t,n||(n=i,i=e.Widget),e.expr[":"][s.toLowerCase()]=function(t){return!!e.data(t,s)},e[u]=e[u]||{},a=e[u][t],o=e[u][t]=function(e,t){return this._createWidget?void(arguments.length&&this._createWidget(e,t)):new o(e,t)},e.extend(o,a,{version:n.version,_proto:e.extend({},n),_childConstructors:[]}),r=new i,r.options=e.widget.extend({},r.options),e.each(n,function(t,n){return e.isFunction(n)?void(l[t]=function(){var e=function(){return i.prototype[t].apply(this,arguments)},s=function(e){return i.prototype[t].apply(this,e)};return function(){var t,i=this._super,a=this._superApply;return this._super=e,this._superApply=s,t=n.apply(this,arguments),this._super=i,this._superApply=a,t}}()):void(l[t]=n)}),o.prototype=e.widget.extend(r,{widgetEventPrefix:a?r.widgetEventPrefix||t:t},l,{constructor:o,namespace:u,widgetName:t,widgetFullName:s}),a?(e.each(a._childConstructors,function(t,i){var n=i.prototype;e.widget(n.namespace+"."+n.widgetName,o,i._proto)}),delete a._childConstructors):i._childConstructors.push(o),e.widget.bridge(t,o),o},e.widget.extend=function(t){for(var n,s,a=i.call(arguments,1),o=0,r=a.length;r>o;o++)for(n in a[o])s=a[o][n],a[o].hasOwnProperty(n)&&void 0!==s&&(e.isPlainObject(s)?t[n]=e.isPlainObject(t[n])?e.widget.extend({},t[n],s):e.widget.extend({},s):t[n]=s);return t},e.widget.bridge=function(t,n){var s=n.prototype.widgetFullName||t;e.fn[t]=function(a){var o="string"==typeof a,r=i.call(arguments,1),l=this;return o?this.each(function(){var i,n=e.data(this,s);return"instance"===a?(l=n,!1):n?e.isFunction(n[a])&&"_"!==a.charAt(0)?(i=n[a].apply(n,r),i!==n&&void 0!==i?(l=i&&i.jquery?l.pushStack(i.get()):i,!1):void 0):e.error("no such method '"+a+"' for "+t+" widget instance"):e.error("cannot call methods on "+t+" prior to initialization; attempted to call method '"+a+"'")}):(r.length&&(a=e.widget.extend.apply(null,[a].concat(r))),this.each(function(){var t=e.data(this,s);t?(t.option(a||{}),t._init&&t._init()):e.data(this,s,new n(a,this))})),l}},e.Widget=function(){},e.Widget._childConstructors=[],e.Widget.prototype={widgetName:"widget",widgetEventPrefix:"",defaultElement:"<div>",options:{disabled:!1,create:null},_createWidget:function(i,n){n=e(n||this.defaultElement||this)[0],this.element=e(n),this.uuid=t++,this.eventNamespace="."+this.widgetName+this.uuid,this.bindings=e(),this.hoverable=e(),this.focusable=e(),n!==this&&(e.data(n,this.widgetFullName,this),this._on(!0,this.element,{remove:function(e){e.target===n&&this.destroy()}}),this.document=e(n.style?n.ownerDocument:n.document||n),this.window=e(this.document[0].defaultView||this.document[0].parentWindow)),this.options=e.widget.extend({},this.options,this._getCreateOptions(),i),this._create(),this._trigger("create",null,this._getCreateEventData()),this._init()},_getCreateOptions:e.noop,_getCreateEventData:e.noop,_create:e.noop,_init:e.noop,destroy:function(){this._destroy(),this.element.unbind(this.eventNamespace).removeData(this.widgetFullName).removeData(e.camelCase(this.widgetFullName)),this.widget().unbind(this.eventNamespace).removeAttr("aria-disabled").removeClass(this.widgetFullName+"-disabled ui-state-disabled"),this.bindings.unbind(this.eventNamespace),this.hoverable.removeClass("ui-state-hover"),this.focusable.removeClass("ui-state-focus")},_destroy:e.noop,widget:function(){return this.element},option:function(t,i){var n,s,a,o=t;if(0===arguments.length)return e.widget.extend({},this.options);if("string"==typeof t)if(o={},n=t.split("."),t=n.shift(),n.length){for(s=o[t]=e.widget.extend({},this.options[t]),a=0;a<n.length-1;a++)s[n[a]]=s[n[a]]||{},s=s[n[a]];if(t=n.pop(),1===arguments.length)return void 0===s[t]?null:s[t];s[t]=i}else{if(1===arguments.length)return void 0===this.options[t]?null:this.options[t];o[t]=i}return this._setOptions(o),this},_setOptions:function(e){var t;for(t in e)this._setOption(t,e[t]);return this},_setOption:function(e,t){return this.options[e]=t,"disabled"===e&&(this.widget().toggleClass(this.widgetFullName+"-disabled",!!t),t&&(this.hoverable.removeClass("ui-state-hover"),this.focusable.removeClass("ui-state-focus"))),this},enable:function(){return this._setOptions({disabled:!1})},disable:function(){return this._setOptions({disabled:!0})},_on:function(t,i,n){var s,a=this;"boolean"!=typeof t&&(n=i,i=t,t=!1),n?(i=s=e(i),this.bindings=this.bindings.add(i)):(n=i,i=this.element,s=this.widget()),e.each(n,function(n,o){function r(){return t||a.options.disabled!==!0&&!e(this).hasClass("ui-state-disabled")?("string"==typeof o?a[o]:o).apply(a,arguments):void 0}"string"!=typeof o&&(r.guid=o.guid=o.guid||r.guid||e.guid++);var l=n.match(/^([\w:-]*)\s*(.*)$/),u=l[1]+a.eventNamespace,h=l[2];h?s.delegate(h,u,r):i.bind(u,r)})},_off:function(t,i){i=(i||"").split(" ").join(this.eventNamespace+" ")+this.eventNamespace,t.unbind(i).undelegate(i),this.bindings=e(this.bindings.not(t).get()),this.focusable=e(this.focusable.not(t).get()),this.hoverable=e(this.hoverable.not(t).get())},_delay:function(e,t){function i(){return("string"==typeof e?n[e]:e).apply(n,arguments)}var n=this;return setTimeout(i,t||0)},_hoverable:function(t){this.hoverable=this.hoverable.add(t),this._on(t,{mouseenter:function(t){e(t.currentTarget).addClass("ui-state-hover")},mouseleave:function(t){e(t.currentTarget).removeClass("ui-state-hover")}})},_focusable:function(t){this.focusable=this.focusable.add(t),this._on(t,{focusin:function(t){e(t.currentTarget).addClass("ui-state-focus")},focusout:function(t){e(t.currentTarget).removeClass("ui-state-focus")}})},_trigger:function(t,i,n){var s,a,o=this.options[t];if(n=n||{},i=e.Event(i),i.type=(t===this.widgetEventPrefix?t:this.widgetEventPrefix+t).toLowerCase(),i.target=this.element[0],a=i.originalEvent)for(s in a)s in i||(i[s]=a[s]);return this.element.trigger(i,n),!(e.isFunction(o)&&o.apply(this.element[0],[i].concat(n))===!1||i.isDefaultPrevented())}},e.each({show:"fadeIn",hide:"fadeOut"},function(t,i){e.Widget.prototype["_"+t]=function(n,s,a){"string"==typeof s&&(s={effect:s});var o,r=s?s===!0||"number"==typeof s?i:s.effect||i:t;s=s||{},"number"==typeof s&&(s={duration:s}),o=!e.isEmptyObject(s),s.complete=a,s.delay&&n.delay(s.delay),o&&e.effects&&e.effects.effect[r]?n[t](s):r!==t&&n[r]?n[r](s.duration,s.easing,a):n.queue(function(i){e(this)[t](),a&&a.call(n[0]),i()})}}),e.widget}),function(e){"function"==typeof define&&define.amd?define("jquery-ui/mouse",["jquery","./widget"],e):e(jQuery)}(function(e){var t=!1;return e(document).mouseup(function(){t=!1}),e.widget("ui.mouse",{version:"1.11.4",options:{cancel:"input,textarea,button,select,option",distance:1,delay:0},_mouseInit:function(){var t=this;this.element.bind("mousedown."+this.widgetName,function(e){return t._mouseDown(e)}).bind("click."+this.widgetName,function(i){return!0===e.data(i.target,t.widgetName+".preventClickEvent")?(e.removeData(i.target,t.widgetName+".preventClickEvent"),i.stopImmediatePropagation(),!1):void 0}),this.started=!1},_mouseDestroy:function(){this.element.unbind("."+this.widgetName),this._mouseMoveDelegate&&this.document.unbind("mousemove."+this.widgetName,this._mouseMoveDelegate).unbind("mouseup."+this.widgetName,this._mouseUpDelegate)},_mouseDown:function(i){if(!t){this._mouseMoved=!1,this._mouseStarted&&this._mouseUp(i),this._mouseDownEvent=i;var n=this,s=1===i.which,a="string"==typeof this.options.cancel&&i.target.nodeName?e(i.target).closest(this.options.cancel).length:!1;return s&&!a&&this._mouseCapture(i)?(this.mouseDelayMet=!this.options.delay,this.mouseDelayMet||(this._mouseDelayTimer=setTimeout(function(){n.mouseDelayMet=!0},this.options.delay)),this._mouseDistanceMet(i)&&this._mouseDelayMet(i)&&(this._mouseStarted=this._mouseStart(i)!==!1,!this._mouseStarted)?(i.preventDefault(),!0):(!0===e.data(i.target,this.widgetName+".preventClickEvent")&&e.removeData(i.target,this.widgetName+".preventClickEvent"),this._mouseMoveDelegate=function(e){return n._mouseMove(e)},this._mouseUpDelegate=function(e){return n._mouseUp(e)},this.document.bind("mousemove."+this.widgetName,this._mouseMoveDelegate).bind("mouseup."+this.widgetName,this._mouseUpDelegate),i.preventDefault(),t=!0,!0)):!0}},_mouseMove:function(t){if(this._mouseMoved){if(e.ui.ie&&(!document.documentMode||document.documentMode<9)&&!t.button)return this._mouseUp(t);if(!t.which)return this._mouseUp(t)}return(t.which||t.button)&&(this._mouseMoved=!0),this._mouseStarted?(this._mouseDrag(t),t.preventDefault()):(this._mouseDistanceMet(t)&&this._mouseDelayMet(t)&&(this._mouseStarted=this._mouseStart(this._mouseDownEvent,t)!==!1,this._mouseStarted?this._mouseDrag(t):this._mouseUp(t)),!this._mouseStarted)},_mouseUp:function(i){return this.document.unbind("mousemove."+this.widgetName,this._mouseMoveDelegate).unbind("mouseup."+this.widgetName,this._mouseUpDelegate),this._mouseStarted&&(this._mouseStarted=!1,i.target===this._mouseDownEvent.target&&e.data(i.target,this.widgetName+".preventClickEvent",!0),this._mouseStop(i)),t=!1,!1},_mouseDistanceMet:function(e){return Math.max(Math.abs(this._mouseDownEvent.pageX-e.pageX),Math.abs(this._mouseDownEvent.pageY-e.pageY))>=this.options.distance},_mouseDelayMet:function(){return this.mouseDelayMet},_mouseStart:function(){},_mouseDrag:function(){},_mouseStop:function(){},_mouseCapture:function(){return!0}})}),function(e){"function"==typeof define&&define.amd?define("jquery-ui/slider",["jquery","./core","./mouse","./widget"],e):e(jQuery)}(function(e){return e.widget("ui.slider",e.ui.mouse,{version:"1.11.4",widgetEventPrefix:"slide",options:{animate:!1,distance:0,max:100,min:0,orientation:"horizontal",range:!1,step:1,value:0,values:null,change:null,slide:null,start:null,stop:null},numPages:5,_create:function(){this._keySliding=!1,this._mouseSliding=!1,this._animateOff=!0,this._handleIndex=null,this._detectOrientation(),this._mouseInit(),this._calculateNewMax(),this.element.addClass("ui-slider ui-slider-"+this.orientation+" ui-widget ui-widget-content ui-corner-all"),this._refresh(),this._setOption("disabled",this.options.disabled),this._animateOff=!1},_refresh:function(){this._createRange(),this._createHandles(),this._setupEvents(),this._refreshValue()},_createHandles:function(){var t,i,n=this.options,s=this.element.find(".ui-slider-handle").addClass("ui-state-default ui-corner-all"),a="<span class='ui-slider-handle ui-state-default ui-corner-all' tabindex='0'></span>",o=[];for(i=n.values&&n.values.length||1,s.length>i&&(s.slice(i).remove(),s=s.slice(0,i)),t=s.length;i>t;t++)o.push(a);this.handles=s.add(e(o.join("")).appendTo(this.element)),this.handle=this.handles.eq(0),this.handles.each(function(t){e(this).data("ui-slider-handle-index",t)})},_createRange:function(){var t=this.options,i="";t.range?(t.range===!0&&(t.values?t.values.length&&2!==t.values.length?t.values=[t.values[0],t.values[0]]:e.isArray(t.values)&&(t.values=t.values.slice(0)):t.values=[this._valueMin(),this._valueMin()]),this.range&&this.range.length?this.range.removeClass("ui-slider-range-min ui-slider-range-max").css({left:"",bottom:""}):(this.range=e("<div></div>").appendTo(this.element),i="ui-slider-range ui-widget-header ui-corner-all"),this.range.addClass(i+("min"===t.range||"max"===t.range?" ui-slider-range-"+t.range:""))):(this.range&&this.range.remove(),this.range=null)},_setupEvents:function(){this._off(this.handles),this._on(this.handles,this._handleEvents),this._hoverable(this.handles),this._focusable(this.handles)},_destroy:function(){this.handles.remove(),this.range&&this.range.remove(),this.element.removeClass("ui-slider ui-slider-horizontal ui-slider-vertical ui-widget ui-widget-content ui-corner-all"),this._mouseDestroy()},_mouseCapture:function(t){var i,n,s,a,o,r,l,u,h=this,d=this.options;return d.disabled?!1:(this.elementSize={width:this.element.outerWidth(),height:this.element.outerHeight()},this.elementOffset=this.element.offset(),i={x:t.pageX,y:t.pageY},n=this._normValueFromMouse(i),s=this._valueMax()-this._valueMin()+1,this.handles.each(function(t){var i=Math.abs(n-h.values(t));(s>i||s===i&&(t===h._lastChangedValue||h.values(t)===d.min))&&(s=i,a=e(this),o=t)}),r=this._start(t,o),r===!1?!1:(this._mouseSliding=!0,this._handleIndex=o,a.addClass("ui-state-active").focus(),l=a.offset(),u=!e(t.target).parents().addBack().is(".ui-slider-handle"),this._clickOffset=u?{left:0,top:0}:{left:t.pageX-l.left-a.width()/2,top:t.pageY-l.top-a.height()/2-(parseInt(a.css("borderTopWidth"),10)||0)-(parseInt(a.css("borderBottomWidth"),10)||0)+(parseInt(a.css("marginTop"),10)||0)},this.handles.hasClass("ui-state-hover")||this._slide(t,o,n),this._animateOff=!0,!0))},_mouseStart:function(){return!0},_mouseDrag:function(e){var t={x:e.pageX,y:e.pageY},i=this._normValueFromMouse(t);return this._slide(e,this._handleIndex,i),!1},_mouseStop:function(e){return this.handles.removeClass("ui-state-active"),this._mouseSliding=!1,this._stop(e,this._handleIndex),this._change(e,this._handleIndex),this._handleIndex=null,this._clickOffset=null,this._animateOff=!1,!1},_detectOrientation:function(){this.orientation="vertical"===this.options.orientation?"vertical":"horizontal"},_normValueFromMouse:function(e){var t,i,n,s,a;return"horizontal"===this.orientation?(t=this.elementSize.width,i=e.x-this.elementOffset.left-(this._clickOffset?this._clickOffset.left:0)):(t=this.elementSize.height,i=e.y-this.elementOffset.top-(this._clickOffset?this._clickOffset.top:0)),n=i/t,n>1&&(n=1),0>n&&(n=0),"vertical"===this.orientation&&(n=1-n),s=this._valueMax()-this._valueMin(),a=this._valueMin()+n*s,this._trimAlignValue(a)},_start:function(e,t){var i={handle:this.handles[t],value:this.value()};return this.options.values&&this.options.values.length&&(i.value=this.values(t),i.values=this.values()),this._trigger("start",e,i)},_slide:function(e,t,i){var n,s,a;this.options.values&&this.options.values.length?(n=this.values(t?0:1),2===this.options.values.length&&this.options.range===!0&&(0===t&&i>n||1===t&&n>i)&&(i=n),i!==this.values(t)&&(s=this.values(),s[t]=i,a=this._trigger("slide",e,{handle:this.handles[t],value:i,values:s}),n=this.values(t?0:1),a!==!1&&this.values(t,i))):i!==this.value()&&(a=this._trigger("slide",e,{handle:this.handles[t],value:i}),a!==!1&&this.value(i))},_stop:function(e,t){var i={handle:this.handles[t],value:this.value()};this.options.values&&this.options.values.length&&(i.value=this.values(t),i.values=this.values()),this._trigger("stop",e,i)},_change:function(e,t){if(!this._keySliding&&!this._mouseSliding){var i={handle:this.handles[t],value:this.value()};this.options.values&&this.options.values.length&&(i.value=this.values(t),i.values=this.values()),this._lastChangedValue=t,this._trigger("change",e,i)}},value:function(e){return arguments.length?(this.options.value=this._trimAlignValue(e),this._refreshValue(),void this._change(null,0)):this._value()},values:function(t,i){var n,s,a;if(arguments.length>1)return this.options.values[t]=this._trimAlignValue(i),this._refreshValue(),void this._change(null,t);if(!arguments.length)return this._values();if(!e.isArray(arguments[0]))return this.options.values&&this.options.values.length?this._values(t):this.value();for(n=this.options.values,s=arguments[0],a=0;a<n.length;a+=1)n[a]=this._trimAlignValue(s[a]),this._change(null,a);this._refreshValue()},_setOption:function(t,i){var n,s=0;switch("range"===t&&this.options.range===!0&&("min"===i?(this.options.value=this._values(0),this.options.values=null):"max"===i&&(this.options.value=this._values(this.options.values.length-1),this.options.values=null)),e.isArray(this.options.values)&&(s=this.options.values.length),"disabled"===t&&this.element.toggleClass("ui-state-disabled",!!i),this._super(t,i),t){case"orientation":this._detectOrientation(),this.element.removeClass("ui-slider-horizontal ui-slider-vertical").addClass("ui-slider-"+this.orientation),this._refreshValue(),this.handles.css("horizontal"===i?"bottom":"left","");break;case"value":this._animateOff=!0,this._refreshValue(),this._change(null,0),this._animateOff=!1;break;case"values":for(this._animateOff=!0,this._refreshValue(),n=0;s>n;n+=1)this._change(null,n);this._animateOff=!1;break;case"step":case"min":case"max":this._animateOff=!0,this._calculateNewMax(),this._refreshValue(),this._animateOff=!1;break;case"range":this._animateOff=!0,this._refresh(),this._animateOff=!1}},_value:function(){var e=this.options.value;return e=this._trimAlignValue(e)},_values:function(e){var t,i,n;if(arguments.length)return t=this.options.values[e],t=this._trimAlignValue(t);if(this.options.values&&this.options.values.length){for(i=this.options.values.slice(),n=0;n<i.length;n+=1)i[n]=this._trimAlignValue(i[n]);return i}return[]},_trimAlignValue:function(e){if(e<=this._valueMin())return this._valueMin();if(e>=this._valueMax())return this._valueMax();var t=this.options.step>0?this.options.step:1,i=(e-this._valueMin())%t,n=e-i;return 2*Math.abs(i)>=t&&(n+=i>0?t:-t),parseFloat(n.toFixed(5))},_calculateNewMax:function(){var e=this.options.max,t=this._valueMin(),i=this.options.step,n=Math.floor(+(e-t).toFixed(this._precision())/i)*i;e=n+t,this.max=parseFloat(e.toFixed(this._precision()))},_precision:function(){var e=this._precisionOf(this.options.step);return null!==this.options.min&&(e=Math.max(e,this._precisionOf(this.options.min))),e},_precisionOf:function(e){var t=e.toString(),i=t.indexOf(".");return-1===i?0:t.length-i-1},_valueMin:function(){return this.options.min},_valueMax:function(){return this.max},_refreshValue:function(){var t,i,n,s,a,o=this.options.range,r=this.options,l=this,u=this._animateOff?!1:r.animate,h={};this.options.values&&this.options.values.length?this.handles.each(function(n){i=(l.values(n)-l._valueMin())/(l._valueMax()-l._valueMin())*100,h["horizontal"===l.orientation?"left":"bottom"]=i+"%",e(this).stop(1,1)[u?"animate":"css"](h,r.animate),l.options.range===!0&&("horizontal"===l.orientation?(0===n&&l.range.stop(1,1)[u?"animate":"css"]({left:i+"%"},r.animate),1===n&&l.range[u?"animate":"css"]({width:i-t+"%"},{queue:!1,duration:r.animate})):(0===n&&l.range.stop(1,1)[u?"animate":"css"]({bottom:i+"%"},r.animate),1===n&&l.range[u?"animate":"css"]({height:i-t+"%"},{queue:!1,duration:r.animate}))),t=i}):(n=this.value(),s=this._valueMin(),a=this._valueMax(),i=a!==s?(n-s)/(a-s)*100:0,
h["horizontal"===this.orientation?"left":"bottom"]=i+"%",this.handle.stop(1,1)[u?"animate":"css"](h,r.animate),"min"===o&&"horizontal"===this.orientation&&this.range.stop(1,1)[u?"animate":"css"]({width:i+"%"},r.animate),"max"===o&&"horizontal"===this.orientation&&this.range[u?"animate":"css"]({width:100-i+"%"},{queue:!1,duration:r.animate}),"min"===o&&"vertical"===this.orientation&&this.range.stop(1,1)[u?"animate":"css"]({height:i+"%"},r.animate),"max"===o&&"vertical"===this.orientation&&this.range[u?"animate":"css"]({height:100-i+"%"},{queue:!1,duration:r.animate}))},_handleEvents:{keydown:function(t){var i,n,s,a,o=e(t.target).data("ui-slider-handle-index");switch(t.keyCode){case e.ui.keyCode.HOME:case e.ui.keyCode.END:case e.ui.keyCode.PAGE_UP:case e.ui.keyCode.PAGE_DOWN:case e.ui.keyCode.UP:case e.ui.keyCode.RIGHT:case e.ui.keyCode.DOWN:case e.ui.keyCode.LEFT:if(t.preventDefault(),!this._keySliding&&(this._keySliding=!0,e(t.target).addClass("ui-state-active"),i=this._start(t,o),i===!1))return}switch(a=this.options.step,n=s=this.options.values&&this.options.values.length?this.values(o):this.value(),t.keyCode){case e.ui.keyCode.HOME:s=this._valueMin();break;case e.ui.keyCode.END:s=this._valueMax();break;case e.ui.keyCode.PAGE_UP:s=this._trimAlignValue(n+(this._valueMax()-this._valueMin())/this.numPages);break;case e.ui.keyCode.PAGE_DOWN:s=this._trimAlignValue(n-(this._valueMax()-this._valueMin())/this.numPages);break;case e.ui.keyCode.UP:case e.ui.keyCode.RIGHT:if(n===this._valueMax())return;s=this._trimAlignValue(n+a);break;case e.ui.keyCode.DOWN:case e.ui.keyCode.LEFT:if(n===this._valueMin())return;s=this._trimAlignValue(n-a)}this._slide(t,o,s)},keyup:function(t){var i=e(t.target).data("ui-slider-handle-index");this._keySliding&&(this._keySliding=!1,this._stop(t,i),this._change(t,i),e(t.target).removeClass("ui-state-active"))}}})}),define("app/user/settings",["jquery","app/ui/layout","jquery-ui/slider","bootstrap"],function(e,t){function i(t){var i=e(this),n=(e(i.attr("data-target")),e(i.attr("root-modal")));return e.ajax({type:i.attr("method"),url:i.attr("action"),data:i.serialize(),success:function(t){e.each(t,function(i,s){var a=e("input[name="+i+"]");a.val(s),a.removeClass("alert alert-danger"),e("#"+i).html(s),n.modal("hide"),t.hasOwnProperty("redirect")&&e(location).attr("href",t.redirect)})},error:function(t){e("#id_errors").empty();var i=JSON.parse(t.responseText);e("input").removeClass("alert alert-danger"),e.each(i,function(t,i){e("input[name="+t+"]").addClass("alert alert-danger"),e("#id_errors").prepend('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>'+i+"</div>")})}}),t.preventDefault(),!1}function n(){e("#id_stream_vector_weight, #id_stream_author_weight, #id_stream_journal_weight, #id_stream_reactivity, #id_trend_doc_weight, #id_trend_altmetric_weight").hide().wrap('<div class="slider"></div>'),e(".slider").each(function(t,i){var n=e(i),s=n.find("input");n.slider({range:!1,animate:!1,min:100*parseFloat(s.data("slider-min")),max:100*parseFloat(s.data("slider-max")),step:100*parseFloat(s.data("slider-step")),value:100*parseFloat(s.val()),slide:function(e,t){var i=t.value/100;n.find("span").tooltip("destroy").tooltip({title:i,animation:!1,container:n.find("span")}).tooltip("show"),s.val(i)},create:function(){n.find("span").tooltip({title:s.val(),animation:!1,container:n.find("span")})}})})}function s(i){var n,s;t.setBusy(),s=e.ajax(i),s.done(function(){n=setInterval(function(){e.getJSON("/user/user-update-step",function(e){if(e.done){if(clearInterval(n),e.hasOwnProperty("redirect"))return void(window.location.href=e.redirect);t.setAvailable()}else t.setBusy("<p><strong>"+e.messages[0]+"</strong></p><p>"+e.messages[1]+"</p>")})},1e3)}),s.fail(function(i){var n=JSON.parse(i.responseText);e.each(n,function(t,i){e("#errors").html(i)}),t.setAvailable()})}e(function(){e("form[data-async]").on("submit",i),n(),e("#update-lib, #update-stream, #update-trend").on("click",function(){s(e(this).attr("action"))})})});