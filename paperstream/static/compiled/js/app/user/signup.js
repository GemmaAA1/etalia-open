define("app/util/utils",["jquery"],function(e){function t(e){var t=null;if(document.cookie&&""!=document.cookie)for(var n=document.cookie.split(";"),i=0;i<n.length;i++){var o=jQuery.trim(n[i]);if(o.substring(0,e.length+1)==e+"="){t=decodeURIComponent(o.substring(e.length+1));break}}return t}function n(e){return/^(GET|HEAD|OPTIONS|TRACE)$/.test(e)}var i=t("csrftoken");return e.ajaxSetup({beforeSend:function(e,t){n(t.type)||this.crossDomain||e.setRequestHeader("X-CSRFToken",i)}}),{toggleClass:function(e,t,n){return"undefined"!=typeof n?(n?e.addClass(t):e.removeClass(t),n):e.hasClass(t)?(e.removeClass(t),!1):(e.addClass(t),!0)},getParameterByName:function(e,t){t||(t=window.location.href),e=e.replace(/[\[\]]/g,"\\$&");var n=new RegExp("[?&]"+e+"(=([^&#]*)|&|#|$)"),i=n.exec(t);return i?i[2]?decodeURIComponent(i[2].replace(/\+/g," ")):"":null},bindLoadingButtons:function(t){t.on("click",".loading-button",function(t){var n=e(t.target).closest(".loading-button");return n.data("busy")?(t.stopPropagation(),void t.stopImmediatePropagation()):(n.data("busy",!0),void n.find(".eai").removeAttr("class").addClass("eai eai-loading"))})},restoreLoadingButton:function(e,t){e.data("busy",!1),e.find(".eai").removeAttr("class").addClass("eai "+t)},popup:function(e,t,n,i){t=t||"etalia-popup",n=n||520,i=i||460;var o=Math.round(window.innerWidth/2-n/2),l=Math.round(window.innerHeight/2-i/2),a="menubar=no,toolbar=no,resizable=yes,scrollbars=yes,width="+n+",height="+i+",top="+l+",left="+o;return console.log(a),window.open(e,t,a)}}}),function(e){"function"==typeof define&&define.amd?define("jquery.mousewheel",["jquery"],e):"object"==typeof exports?module.exports=e:e(jQuery)}(function(e){function t(t){var a=t||window.event,s=r.call(arguments,1),u=0,d=0,h=0,f=0,p=0,g=0;if(t=e.event.fix(a),t.type="mousewheel","detail"in a&&(h=-1*a.detail),"wheelDelta"in a&&(h=a.wheelDelta),"wheelDeltaY"in a&&(h=a.wheelDeltaY),"wheelDeltaX"in a&&(d=-1*a.wheelDeltaX),"axis"in a&&a.axis===a.HORIZONTAL_AXIS&&(d=-1*h,h=0),u=0===h?d:h,"deltaY"in a&&(h=-1*a.deltaY,u=h),"deltaX"in a&&(d=a.deltaX,0===h&&(u=-1*d)),0!==h||0!==d){if(1===a.deltaMode){var v=e.data(this,"mousewheel-line-height");u*=v,h*=v,d*=v}else if(2===a.deltaMode){var b=e.data(this,"mousewheel-page-height");u*=b,h*=b,d*=b}if(f=Math.max(Math.abs(h),Math.abs(d)),(!l||l>f)&&(l=f,i(a,f)&&(l/=40)),i(a,f)&&(u/=40,d/=40,h/=40),u=Math[u>=1?"floor":"ceil"](u/l),d=Math[d>=1?"floor":"ceil"](d/l),h=Math[h>=1?"floor":"ceil"](h/l),c.settings.normalizeOffset&&this.getBoundingClientRect){var m=this.getBoundingClientRect();p=t.clientX-m.left,g=t.clientY-m.top}return t.deltaX=d,t.deltaY=h,t.deltaFactor=l,t.offsetX=p,t.offsetY=g,t.deltaMode=0,s.unshift(t,u,d,h),o&&clearTimeout(o),o=setTimeout(n,200),(e.event.dispatch||e.event.handle).apply(this,s)}}function n(){l=null}function i(e,t){return c.settings.adjustOldDeltas&&"mousewheel"===e.type&&t%120===0}var o,l,a=["wheel","mousewheel","DOMMouseScroll","MozMousePixelScroll"],s="onwheel"in document||document.documentMode>=9?["wheel"]:["mousewheel","DomMouseScroll","MozMousePixelScroll"],r=Array.prototype.slice;if(e.event.fixHooks)for(var u=a.length;u;)e.event.fixHooks[a[--u]]=e.event.mouseHooks;var c=e.event.special.mousewheel={version:"3.1.12",setup:function(){if(this.addEventListener)for(var n=s.length;n;)this.addEventListener(s[--n],t,!1);else this.onmousewheel=t;e.data(this,"mousewheel-line-height",c.getLineHeight(this)),e.data(this,"mousewheel-page-height",c.getPageHeight(this))},teardown:function(){if(this.removeEventListener)for(var n=s.length;n;)this.removeEventListener(s[--n],t,!1);else this.onmousewheel=null;e.removeData(this,"mousewheel-line-height"),e.removeData(this,"mousewheel-page-height")},getLineHeight:function(t){var n=e(t),i=n["offsetParent"in e.fn?"offsetParent":"parent"]();return i.length||(i=e("body")),parseInt(i.css("fontSize"),10)||parseInt(n.css("fontSize"),10)||16},getPageHeight:function(t){return e(t).height()},settings:{adjustOldDeltas:!0,normalizeOffset:!0}};e.fn.extend({mousewheel:function(e){return e?this.bind("mousewheel",e):this.trigger("mousewheel")},unmousewheel:function(e){return this.unbind("mousewheel",e)}})}),define("app/ui/layout/flap",["jquery","jquery.mousewheel"],function(e){var t=function(t){t=t||{},this.config=e.extend({debug:!0,flap:null,side:null,button:null,backdrop:null,mobileMaxWidth:992},t),this.$flap=e(this.config.flap),this.$button=e(this.config.button),this.$backdrop=e(this.config.backdrop),this.mobile=null,this.clickHandlersEnabled=!1,this.scrollHandlersEnabled=!1,this.opened=!1,this.scrollTop=0;var n=this;this.log=function(e){n.config.debug&&console.log("[Flap "+this.config.side+"] "+e)},this.buttonClickHandler=function(){n.log("Button click handler"),n.open()},this.backdropClickHandler=function(){n.log("Overlay click handler"),n.close()},this.flapMouseWheelHandler=function(e){return n.log("Window scroll handler"),n.affix(50*e.deltaY),e.preventDefault(),!1},this.backdropMouseWheelHandler=function(e){return e.preventDefault(),!1}};return t.prototype.init=function(){this.log("Flap init");var e=window.innerWidth<this.config.mobileMaxWidth;e!=this.mobile&&(this.close(),this.disableScrollHandlers(),e?this.enableClickHandlers():(this.disableClickHandlers(),this.enableScrollHandlers())),this.mobile=e},t.prototype.enableClickHandlers=function(){this.clickHandlersEnabled||(this.log("enableClickHandlers"),this.$button.on("click",this.buttonClickHandler),this.$backdrop.on("click",this.backdropClickHandler),this.clickHandlersEnabled=!0)},t.prototype.disableClickHandlers=function(){this.clickHandlersEnabled&&(this.log("disableClickHandlers"),this.$button.off("click",this.buttonClickHandler),this.$backdrop.off("click",this.backdropClickHandler),this.clickHandlersEnabled=!1)},t.prototype.enableScrollHandlers=function(){this.scrollHandlersEnabled||(this.log("enableScrollHandlers"),this.$flap.on("mousewheel",this.flapMouseWheelHandler),this.$backdrop.on("mousewheel",this.flapMouseWheelHandler),this.scrollHandlersEnabled=!0)},t.prototype.disableScrollHandlers=function(){this.scrollHandlersEnabled&&(this.log("disableScrollHandlers"),this.$flap.off("mousewheel",this.flapMouseWheelHandler),this.$backdrop.off("mousewheel",this.flapMouseWheelHandler),this.scrollHandlersEnabled=!1)},t.prototype.open=function(){return this.opened?this:(this.opened=!0,e("body").addClass("flap-opened").addClass(this.config.side),this.$flap.css({top:0}),void this.enableScrollHandlers())},t.prototype.close=function(){return this.opened?(e("body").removeClass("flap-opened").removeClass(this.config.side),this.$flap.css({top:0}),this.disableScrollHandlers(),void(this.opened=!1)):this},t.prototype.affix=function(e){var t=parseInt(this.$flap.css("top"))+e,n=window.innerHeight-this.$flap.outerHeight();t>0?t=0:n>t&&(t=n),this.$flap.css({top:t}),this.scrollTop=window.scrollY},t}),define("app/ui/layout/interact",["jquery"],function(e){var t=function(t){this.config=e.extend({debug:!1},t)};return t.prototype.log=function(){return this.config.debug&&console.log("[Interact] "+arguments[0],Array.prototype.splice.call(arguments,1)),this},t.prototype.init=function(){this.log("init()");var e=this,t=0,n=setInterval(function(){return e.log("check ga",t,typeof ga),"function"==typeof ga?(e.initHandlers(),void clearInterval(n)):(t++,void(t>20&&clearInterval(n)))},500);return this},t.prototype.initHandlers=function(){this.log("initHandlers()");var t=this,n=e("body"),i="control";return n.on("etalia.control.search.change",function(e,n){var o={hitType:"event",eventCategory:i,eventAction:"search",eventLabel:n.expression};t.log("ga.send()",o),ga("send",o)}).on("etalia.control.timespan.change",function(e,n){var o={hitType:"event",eventCategory:i,eventAction:"timespan",eventLabel:n.label,eventValue:n.value};t.log("ga.send()",o),ga("send",o)}).on("etalia.control.cluster.change",function(e,n){var o={hitType:"event",eventCategory:i,eventAction:"cluster",eventLabel:n.label,eventValue:n.value};t.log("ga.send()",o),ga("send",o)}).on("etalia.control.pinned.change",function(e,n){var o={hitType:"event",eventCategory:i,eventAction:n.active?"pin":"unpin"};t.log("ga.send()",o),ga("send",o)}).on("etalia.control.filters.change",function(e,n){var o={hitType:"event",eventCategory:i,eventAction:(n.active?"select_":"deselect_")+n.group,eventLabel:n.label,eventValue:n.value};t.log("ga.send()",o),ga("send",o)}),this},(new t).init()}),define("app/ui/layout/invite",["jquery"],function(e){function t(t){t.preventDefault();var i=n.find(".modal-form"),o=n.find(".modal-result");return e.ajax({type:i.attr("method"),url:i.attr("action"),data:i.serialize(),success:function(e){i.hide().find("input, textarea").val(""),o.show()},error:function(t){console.log(t.responseText);var i=n.find(".form-errors").empty(),o=JSON.parse(t.responseText);e.each(o,function(e,t){i.prepend('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>'+t+"</div>")})}}),!1}var n;e(function(){n=e("#invite-modal"),e("#open-invite").on("click",function(){e("#profile-dropdown").hide(),n.modal("show")}),n.find("form[data-async]").on("submit",t)})}),define("app/ui/layout/close-alerts",["jquery"],function(e){e(function(){e("a.close[close-href]").click(function(t){t.preventDefault(),e.post(e(this).attr("close-href"),"",function(){})})})}),define("app/ui/layout",["jquery","app/util/utils","app/ui/layout/flap","app/ui/layout/interact","app/ui/layout/invite","app/ui/layout/close-alerts"],function(e,t,n){var i=function(t){this.config=e.extend({debug:!1,leftFlap:"#nav-flap",leftFlapButton:"#toggle-nav",rightFlap:"#filter-flap",rightFlapButton:"#toggle-filter",backdrop:"#backdrop"},t),this.leftFlap=null,this.rightFlap=null,this.busy=!1;var n=this;this.log=function(e){n.config.debug&&console.log("[Layout] "+e)},this.resizeTimeout=null,this.resizeHandler=function(){n.log("Resize handler"),n.leftFlap&&n.leftFlap.init(),n.rightFlap&&n.rightFlap.init()}};i.prototype.init=function(){var i=e(this.config.leftFlap),o=e(this.config.leftFlapButton),l=e(this.config.rightFlap),a=e(this.config.rightFlapButton),s=!1;if(this.log("Init"),i.length?(this.log("Left flap found"),this.leftFlap=new n({debug:this.config.debug,flap:i,side:"left",button:o,backdrop:this.config.backdrop,mobileMaxWidth:992}),s=!0):o.hide(),l.length?(this.log("Right flap found"),this.rightFlap=new n({debug:this.config.debug,flap:l,side:"right",button:this.config.rightFlapButton,backdrop:this.config.backdrop,mobileMaxWidth:1200}),s=!0):a.hide(),s){var r=this;e(window).on("resize",function(){r.resizeTimeout&&clearTimeout(r.resizeTimeout),r.resizeTimeout=setTimeout(r.resizeHandler,100)}),r.resizeHandler()}var u=e("#toggle-profile"),c=e("#profile-dropdown");c.length&&(u.on("click",function(n){t.toggleClass(e(n.delegateTarget),"active")?c.show():c.hide()}),c.on("click",function(e){e.stopPropagation()}),e(window).on("click",function(t){0==e(t.target).closest("#toggle-profile").length&&(e("#toggle-profile").removeClass("active"),e("#profile-dropdown").hide())}))},i.prototype.setBusy=function(t){e("#busy-content").html(t),this.busy||(this.busy=!0,e("body").addClass("busy"))},i.prototype.setAvailable=function(){this.busy&&(this.busy=!1,e("#busy-content").html(""),e("body").removeClass("busy"))};var o=new i;return o.init(),o}),define("app/user/signup",["jquery","app/ui/layout","bootstrap"],function(e,t){function n(n){n.preventDefault(),t.setBusy();var i=e(this),o=e.ajax({url:i.attr("action"),type:i.attr("method"),data:i.serialize(),dataType:"json"});return o.done(function(n){e("#id_errors").empty(),n.hasOwnProperty("redirect")?window.location.href=n.redirect:t.setAvailable()}),o.fail(function(n){e("#id_errors").empty(),i.find("input").removeClass("alert alert-danger");var o=JSON.parse(n.responseText);e.each(o,function(e,t){i.find("input[name="+e+"]").addClass("alert alert-danger"),i.find("#id_errors").prepend('<div class="error-message"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>'+t+"</div>")}),t.setAvailable()}),!1}e(function(){e("form[data-async]").on("submit",n)})});