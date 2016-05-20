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

        this.$nextButton = $('.detail-nav-left button');
        this.$prevButton = $('.detail-nav-right button');

        this.loaded = false;
        this.loadXhr = null;

        this.paper = new Paper(options).init();

        var that = this;
        $('.detail-nav-center button, #backdrop').on('click', function() {
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
                if (that.paper.getId() == result.getId() && result.isBanned()) {
                    that.close();
                }
            })
            .on('etalia.detail.loaded', function() {
                that.paper.update();
            });

        // Close on click out
        this.$document
            .on('click', function(e) {
                if (0 < $(e.target).closest('.inner').length) {
                    return;
                }
                that.close();
            })
            .on('click', '.neighbors-thumbs .thumb .title a', function(e) {
                that.load($(e.target).closest('.thumb'));

                e.preventDefault();
                return false;
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
            var $prevTitle = $prev.find('.title').clone();
            $prevTitle.find('span').remove();

            that.$prevButton.show().attr('title', $prevTitle.text());
            that.$prevButton.find('> span').html($prevTitle.text());
            that.$prevButton.find('> button').on('click', function() {
                that.load($prev);
            });
        }
        // Next button
        if ($next.length) {
            var $nextTitle = $next.find('.title').clone();
            $nextTitle.find('span').remove();

            that.$nextButton.show().attr('title', $nextTitle.text());
            that.$nextButton.find('> span').html($nextTitle.text());
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
