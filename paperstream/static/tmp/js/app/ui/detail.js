define(['jquery', 'app/util/sticky'], function ($, Sticky) {

    var Detail = function (options) {

        this.config = $.extend({
            debug: false,
            element: '#detail',
            document: '.document',
            actions: '.inner > .actions'
        }, options);

        this.$element = $(this.config.element);

        this.$document = this.$element.find(this.config.document);
        this.actions = null;

        this.$nextButton = $('#detail-next');
        this.$prevButton = $('#detail-prev');

        this.loaded = false;
        this.loadXhr = null;

        var that = this;
        $('#detail-close, #backdrop').on('click', function() {
            that.close();
        });
    };
    Detail.prototype.load = function ($thumb) {
        var that = this,
            $prev = $thumb.prev(),
            $next = $thumb.next();

        $(that).trigger('etalia.detail.loading');

        that.clear();

        that.$document.hide();
        $('body').addClass('detail-opened');

        // Previous button
        if ($prev.length) {
            var prevTitle = $prev.find('.title').text();
            that.$prevButton.show().attr('title', prevTitle);
            that.$prevButton.find('> span').html(prevTitle);
            that.$prevButton.find('> button').on('click', function() {
                that.load($prev);
            });
        }
        // Next button
        if ($next.length) {
            var nextTitle = $next.find('.title').text();
            that.$nextButton.show().attr('title', nextTitle);
            that.$nextButton.find('> span').html(nextTitle);
            that.$nextButton.find('> button').on('click', function() {
                that.load($next);
            });
        }

        var uri = $thumb.find('.title a').attr('href');
        $.get(uri, function(html) {
            that.$document.find('.inner').html(html);

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

            $(that).trigger('etalia.detail.loaded');
        });
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

        var $navButtons = $('.detail-nav');
        $navButtons.hide().removeAttr('title');
        $navButtons.find('> button').off('click');
        $navButtons.find('> span').empty();
    };
    Detail.prototype.close = function () {
        this.clear();

        $('body').removeClass('detail-opened');
    };

    return Detail;
});
