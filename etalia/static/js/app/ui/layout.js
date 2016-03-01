define([
    'jquery',
    'app/util/utils',
    'app/ui/layout/flap',
    'app/ui/layout/interact',
    'app/ui/layout/invite',
    'app/ui/layout/close-alerts'
], function ($, Util, Flap) {

    var Layout = function (config) {
        this.config = $.extend({
            debug: false,
            leftFlap: '#nav-flap',
            leftFlapButton: '#toggle-nav',
            rightFlap: '#filter-flap',
            rightFlapButton: '#toggle-filter',
            backdrop: '#backdrop'
        }, config);

        this.leftFlap = null;
        this.rightFlap = null;

        this.busy = false;

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
            $leftFlapButton = $(this.config.leftFlapButton),
            $rightFlap = $(this.config.rightFlap),
            $rightFlapButton = $(this.config.rightFlapButton),
            attacheResizeHandler = false;

        this.log('Init');

        // Left Flap
        if ($leftFlap.length) {
            this.log('Left flap found');
            this.leftFlap = new Flap({
                debug: this.config.debug,
                flap: this.config.leftFlap,
                side: 'left',
                button: this.config.leftFlapButton,
                backdrop: this.config.backdrop,
                mobileMaxWidth: 992
            });
            attacheResizeHandler = true;
        } else {
            $leftFlapButton.hide();
        }

        // Right flap
        if ($rightFlap.length) {
            this.log('Right flap found');
            this.rightFlap = new Flap({
                debug: this.config.debug,
                flap: this.config.rightFlap,
                side: 'right',
                button: this.config.rightFlapButton,
                backdrop: this.config.backdrop,
                mobileMaxWidth: 1200
            });
            attacheResizeHandler = true;
        } else {
            $rightFlapButton.hide();
        }

        // Resize handler
        if (attacheResizeHandler) {
            var that = this;
            $(window).on('resize', function () {
                if (that.resizeTimeout) {
                    clearTimeout(that.resizeTimeout);
                }
                that.resizeTimeout = setTimeout(that.resizeHandler, 100);
            });
            that.resizeHandler();
        }

        // Profile dropdown
        var $toggleProfile = $('#toggle-profile'),
            $profileDropDown = $('#profile-dropdown');

        if ($profileDropDown.length) {
            $toggleProfile.on('click', function (e) {
                if (Util.toggleClass($(e.delegateTarget), 'active')) {
                    $profileDropDown.show();
                } else {
                    $profileDropDown.hide();
                }
            });
            $profileDropDown.on('click', function (e) {
                e.stopPropagation();
            });

            $(window).on('click', function(e) {
                // Profile dropdown
                if (0 == $(e.target).closest('#toggle-profile').length) {
                    $('#toggle-profile').removeClass('active');
                    $('#profile-dropdown').hide();
                }
            });
        }
    };

    Layout.prototype.setBusy = function(content) {
        $('#busy-content').html(content);

        if (this.busy) {
            return;
        }
        this.busy = true;

        $('body').addClass('busy');
    };

    Layout.prototype.setAvailable = function() {
        if (!this.busy) {
            return;
        }
        this.busy = false;

        $('#busy-content').html('');
        $('body').removeClass('busy');
    };

    var layout = new Layout();
    layout.init();

    return layout;
});
