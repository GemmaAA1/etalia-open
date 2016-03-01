define(['jquery'], function($) {

    var Sticky = function(options) {
        var config = $.extend({
            debug: true,
            element: null,
            parent: null,
            top: 0,
            bottom: 0
        }, options);

        // TODO check config

        var $element = $(config.element),
            $parent = config.parent ? $(config.parent) : $(window),
            resizeTimeout = null,
            lastScrollTop = null,
            elementHeight = 0,
            parentHeight = 0,
            parentScrollTop = 0,
            elementTop = 0,
            overflow = false;

        function log(msg) {
            if (config.debug) {
                console.log('[Sticky] ' + msg);
            }
        }

        // Scroll handler
        function scrollHandler() {
            parentScrollTop = $parent.scrollTop();
            elementTop = 0;

            log('scrollHandler : ' + elementTop + ' / ' + parentScrollTop);

            if (null === lastScrollTop || !overflow) {
                if (elementTop < parentScrollTop + config.top) {
                    elementTop = parentScrollTop + config.top;
                }
                log('fixed : ' + elementTop);
            } else {
                elementTop = parseInt($element.css('top'));
                var scrollDelta = lastScrollTop - parentScrollTop;
                elementTop += scrollDelta;

                if (0 > scrollDelta) {  // Scrolling down
                    if (elementTop + elementHeight < parentScrollTop + parentHeight - config.bottom) {
                        elementTop = parentScrollTop + parentHeight - elementHeight - config.bottom;
                    }
                    log('scrolling down : ' + elementTop);
                } else {                // Scrolling up
                    if (elementTop > parentScrollTop + config.top) {
                        elementTop = parentScrollTop + config.top;
                    }
                    log('scrolling up : ' + elementTop);
                }
            }

            lastScrollTop = parentScrollTop;
            $element.css({'top': elementTop});
        }

        // Resize handler
        function resizeTimeoutHandler() {
            elementHeight = $element.outerHeight();
            parentHeight = $parent.innerHeight();

            log('resizeTimeoutHandler : ' + elementHeight + ' / ' + parentHeight);

            overflow = elementHeight > parentHeight;
        }
        function resizeHandler() {
            if (resizeTimeout) {
                clearTimeout(resizeTimeout);
            }
            resizeTimeout = setTimeout(resizeTimeoutHandler, 150);
        }

        this.addHandlers = function() {
            $parent.on('scroll', scrollHandler);
            $(window).on('resize', resizeHandler);
            resizeTimeoutHandler();
            scrollHandler();
        };
        this.removeHandlers = function() {
            $parent.off('scroll', scrollHandler);
            $(window).off('resize', resizeHandler);
            $element.css({position: 'absolute', top: '0px'});
            lastScrollTop = null;
        };

        this.enabled = false;
    };

    Sticky.prototype.enable = function() {
        if (this.enabled) return;
        this.addHandlers();
        this.enabled = true;
    };
    Sticky.prototype.disable = function() {
        if (!this.enabled) return;
        this.removeHandlers();
        this.enabled = false;
    };

    return Sticky;
});
