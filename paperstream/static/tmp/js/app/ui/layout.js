define(['jquery', 'app/ui/flap'], function ($, Flap) {

    var Layout = function (config) {
        this.config = $.extend({
            debug: true,
            leftFlap: '#nav-flap',
            leftFlapButton: '#toggle-nav',
            rightFlap: '#filter-flap',
            rightFlapButton: '#toggle-filter',
            backdrop: '#backdrop'
        }, config);

        this.leftFlap = null;
        this.rightFlap = null;

        var that = this;
        this.log = function(msg) {
            if (that.config.debug) {
                console.log('[Layout] ' + msg);
            }
        };

        this.resizeTimeout = null;
        this.resizeHandler = function () {
            that.log('Resize handler');
            if (that.leftFlap) {
                that.leftFlap.init();
            }
            if (that.rightFlap) {
                that.rightFlap.init();
            }
        };
    };
    Layout.prototype.init = function () {
        var $leftFlap = $(this.config.leftFlap),
            $rightFlap = $(this.config.rightFlap);

        this.log('Init');

        if ($leftFlap.length) {
            this.log('Left flap found');
            this.leftFlap = new Flap({
                debug: this.config.debug,
                flap: $leftFlap,
                side: 'left',
                button: this.config.leftFlapButton,
                backdrop: this.config.backdrop,
                mobileMaxWidth: 992
            });
        }
        if ($rightFlap.length) {
            this.log('Right flap found');
            this.rightFlap = new Flap({
                debug: this.config.debug,
                flap: $rightFlap,
                side: 'right',
                button: this.config.rightFlapButton,
                backdrop: this.config.backdrop,
                mobileMaxWidth: 1200
            });
        }

        var that = this;
        $(window).on('resize', function () {
            if (that.resizeTimeout) {
                clearTimeout(that.resizeTimeout);
            }
            that.resizeTimeout = setTimeout(that.resizeHandler, 100);
        });
        that.resizeHandler();
    };

    return Layout;
});
