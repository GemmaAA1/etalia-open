define(['jquery'], function($) {

    var Flap = function(options) {
        options = options || {};

        this.config = $.extend({
            debug: true,
            flap: null,
            side: null,
            button: null,
            backdrop: null,
            mobileMaxWidth: 992
        }, options);

        this.$flap = $(this.config.flap);
        this.$button = $(this.config.button);
        this.$backdrop = $(this.config.backdrop);

        this.mobile = null;
        this.clickHandlersEnabled = false;
        this.scrollHandlersEnabled = false;
        this.opened = false;
        this.scrollTop = 0;

        var that = this;
        this.log = function(msg) {
            if (that.config.debug) {
                console.log('[Flap ' + this.config.side + '] ' + msg);
            }
        };
        this.buttonClickHandler = function() {
            that.log('Button click handler');
            that.open();
        };
        this.backdropClickHandler = function() {
            that.log('Overlay click handler');
            that.close();
        };
        this.windowScrollHandler = function() {
            that.log('Window scroll handler');
            that.affix();
        };
    };
    Flap.prototype.init = function() {
        var mobile = window.innerWidth < this.config.mobileMaxWidth;
        if (mobile != this.mobile) {
            this.close();
            this.disableScrollHandlers();
            if (mobile) {
                this.enableClickHandlers();
            } else {
                this.disableClickHandlers();
                this.enableScrollHandlers();
            }
        }
        this.mobile = mobile;
    };
    Flap.prototype.enableClickHandlers = function() {
        if (this.clickHandlersEnabled) {
            return;
        }

        this.$button.on('click', this.buttonClickHandler);
        this.$backdrop.on('click', this.backdropClickHandler);

        this.clickHandlersEnabled = true;
    };
    Flap.prototype.disableClickHandlers = function() {
        if (!this.clickHandlersEnabled) {
            return;
        }

        this.$button.off('click', this.buttonClickHandler);
        this.$backdrop.off('click', this.backdropClickHandler);

        this.clickHandlersEnabled = false;
    };
    Flap.prototype.enableScrollHandlers = function() {
        if (this.scrollHandlersEnabled) {
            return;
        }

        $(window).on('scroll', this.windowScrollHandler);

        this.scrollHandlersEnabled = true;
    };
    Flap.prototype.disableScrollHandlers = function() {
        if (!this.scrollHandlersEnabled) {
            return;
        }

        $(window).off('scroll', this.windowScrollHandler);

        this.scrollHandlersEnabled = false;
    };
    Flap.prototype.open = function() {
        if (this.opened) {
            return this;
        }

        this.opened = true;

        $('body')
            .addClass('flap-opened')
            .addClass(this.config.side);

        this.$flap.css({'top': 0});
        this.enableScrollHandlers();
    };
    Flap.prototype.close = function() {
        if (!this.opened) {
            return this;
        }

        $('body')
            .removeClass('flap-opened')
            .removeClass(this.config.side);

        this.$flap.css({'top': 0});
        this.disableScrollHandlers();

        this.opened = false;
    };
    Flap.prototype.affix = function() {
        var top = parseInt(this.$flap.css('top')) + this.scrollTop - window.scrollY,
            min = window.innerHeight - this.$flap.outerHeight();

        if (top > 0) {
            top = 0;
        } else if (top < min) {
            top = min;
        }

        // TODO Fix flap top when to high in page

        this.$flap.css({'top': top});
        this.scrollTop = window.scrollY;
    };

    return Flap;
});
