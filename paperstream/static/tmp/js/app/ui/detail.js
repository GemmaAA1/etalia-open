define(['jquery', 'app/util/sticky'], function ($, Sticky) {

    var Detail = function (options) {

        this.config = $.extend({
            debug: false,
            element: '#detail',
            document: '.document',
            loading: '.loading',
            actions: '.inner > .actions'
        }, options);

        this.$element = $(this.config.element);

        this.$document = this.$element.find(this.config.document);
        this.$loading = this.$element.find(this.config.loading);
        this.actions = null;

        this.loaded = false;
        this.loadXhr = null;
    };
    Detail.prototype.load = function () {
        var that = this;
        that.clear();

        that.$document.hide();
        that.$loading.show();

        $('body').addClass('detail-opened');

        // Fake XHR
        that.loadXhr = setTimeout(function () {
            that.$loading.hide();
            that.$document.show();

            that.actions = new Sticky({
                debug: that.config.debug,
                element: that.$element.find(that.config.actions),
                parent: that.$document,
                top: 20,
                bottom: 20
            });
            that.actions.enable();

            that.loaded = true;
            that.loadXhr = null;
        }, 1500);
    };
    Detail.prototype.clear = function () {
        this.loaded = false;
        if (this.loadXhr) {
            clearTimeout(this.loadXhr);
            this.loadXhr = null;
        }
        if (this.actions) {
            this.actions.disable();
            this.actions = null;
        }
    };
    Detail.prototype.close = function () {
        this.clear();

        $('body').removeClass('detail-opened');
    };

    return Detail;
});
