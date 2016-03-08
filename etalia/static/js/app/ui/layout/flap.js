define(['jquery'], function($) {

    var Flap = function(options) {
        options = options || {};

        this.config = $.extend({
            debug: false,
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
        this.opened = false;

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

        this.$flap.on('redraw', function() {
            that.init();
        });
    };
    Flap.prototype.init = function() {
        this.log('Flap init');
        this.$flap.css({height: window.innerHeight + 'px'});
        var mobile = window.innerWidth < this.config.mobileMaxWidth;
        if (mobile != this.mobile) {
            this.close();
            if (mobile) {
                this.enableClickHandlers();
            } else {
                this.disableClickHandlers();
            }
        }
        this.mobile = mobile;
    };
    Flap.prototype.enableClickHandlers = function() {
        if (this.clickHandlersEnabled) {
            return;
        }
        this.log('enableClickHandlers');

        this.$button.on('click', this.buttonClickHandler);
        this.$backdrop.on('click', this.backdropClickHandler);

        this.clickHandlersEnabled = true;
    };
    Flap.prototype.disableClickHandlers = function() {
        if (!this.clickHandlersEnabled) {
            return;
        }
        this.log('disableClickHandlers');

        this.$button.off('click', this.buttonClickHandler);
        this.$backdrop.off('click', this.backdropClickHandler);

        this.clickHandlersEnabled = false;
    };
    Flap.prototype.open = function() {
        if (this.opened) {
            return this;
        }

        $('body')
            .addClass('flap-opened')
            .addClass(this.config.side);

        this.opened = true;
    };
    Flap.prototype.close = function() {
        if (!this.opened) {
            return this;
        }

        $('body')
            .removeClass('flap-opened')
            .removeClass(this.config.side);

        this.opened = false;
    };

    return Flap;
});
