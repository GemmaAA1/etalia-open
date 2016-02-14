define(['jquery', 'app/ui/paper'], function ($, Paper) {

    var Detail = function (options) {

        this.config = $.extend({
            debug: false,
            element: '#detail',
            document: '.document',
            actions: '.inner > .actions'
        }, options);

        this.$element = $(this.config.element);
        this.$document = this.$element.find(this.config.document);

        this.$nextButton = $('#detail-next');
        this.$prevButton = $('#detail-prev');

        this.loaded = false;
        this.loadXhr = null;

        this.paper = new Paper(options).init();

        var that = this;
        $('#detail-close, #backdrop').on('click', function() {
            that.close();
        });
    };

    Detail.prototype.log = function() {
        if (this.config.debug) {
            console.log('[Detail] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }

        return this;
    };

    Detail.prototype.init = function() {
        var that = this;

        // API events
        $('body')
            .on('etalia.publication.ban', function(e, result) {
                if (that.id == result.getId() && result.isBanned()) {
                    that.close();
                }
            })
            .on('etalia.detail.loaded', function() {
                that.paper.update();
            });

        return this;
    };

    Detail.prototype.load = function ($thumb) {
        var that = this,
            $body = $('body'),
            $prev = $thumb.prev(),
            $next = $thumb.next();

        that.clear();

        $body
            .addClass('detail-opened')
            .trigger('etalia.detail.loading');

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

        // Paper async loading
        this.loadXhr = $.ajax({
            method: 'GET',
            url:    $thumb.find('.title a').attr('href')
        })
        .done(function(html) {
            that.$document.find('.inner').html(html);

            that.$document.show();

            that.loaded = true;
            that.loadXhr = null;

            $body.trigger('etalia.detail.loaded');
        })
        .fail(function(xrh, status, error) {
            that.log('Load failure', xrh, status, error);
        });

        return this;
    };

    Detail.prototype.clear = function () {
        this.loaded = false;
        if (this.loadXhr) {
            this.loadXhr.abort();
            this.loadXhr = null;
        }

        this.paper.clear();

        this.$document.hide().find('.inner').html('');

        var $navButtons = $('.detail-nav');
        $navButtons.hide().removeAttr('title');
        $navButtons.find('> button').off('click');
        $navButtons.find('> span').empty();

        return this;
    };

    Detail.prototype.close = function () {
        this.clear();

        $('body').removeClass('detail-opened');

        return this;
    };

    return Detail;
});
