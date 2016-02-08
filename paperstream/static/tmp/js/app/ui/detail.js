define(['jquery', 'app/api', 'app/util/sticky'], function ($, Api, Sticky) {

    var Detail = function (options) {

        this.config = $.extend({
            debug: false,
            element: '#detail',
            document: '.document',
            actions: '.inner > .actions'
        }, options);

        this.$element = $(this.config.element);

        this.$document = this.$element.find(this.config.document);
        this.paperId = null;
        this.$actions = null;
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

    Detail.prototype.init = function() {
        var that = this;

        // Pin button
        this.$document.on('click', '.detail-pin', function(e) {
            if (!that.paperId) throw 'Undefined paper id';

            var $button = $(e.target).closest('.detail-pin');

            Api.pin(that.paperId) // TODO source = window.location.pathname ?
                .done(function(pinned) {
                    $button.toggleClass('active', pinned);
                })
                .fail(function(error) {
                    console.log(error);
                });

            e.stopPropagation();
            return false;
        });

        // Ban button
        this.$document.on('click', '.detail-ban', function(e) {
            if (!that.paperId) throw 'Undefined paper id';

            Api.ban(that.paperId) // TODO source = window.location.pathname ?
                .done(function(banned) {
                    if (banned) {
                        that.close();
                    }
                })
                .fail(function(error) {
                    console.log(error);
                });

            e.stopPropagation();
            return false;
        });

        // Add to library button
        this.$document.on('click', '.detail-library-add', function(e) {
            if (!that.paperId) throw 'Undefined paper id';

            var $button = $(e.target).closest('.detail-library-add');
            // TODO Don't add twice ?
            if ($button.hasClass('active')) {
                return;
            }

            Api.add(that.paperId)
                .done(function(added) {
                    $button.toggleClass('active', added);
                })
                .fail(function(error) {
                    console.log(error);
                });

            e.stopPropagation();
            return false;
        });

        return this;
    };

    Detail.prototype.load = function ($thumb) {
        var that = this,
            $prev = $thumb.prev(),
            $next = $thumb.next();

        $(that).trigger('etalia.detail.loading');

        that.clear();

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

        // Paper async loading
        var uri = $thumb.find('.title a').attr('href');
        this.loadXhr = $.get(uri)
            .done(function(html) {
                that.$document.find('.inner').html(html);

                that.$actions = that.$element.find(that.config.actions);
                that.paperId = parseInt(that.$actions.data('paper-id'));

                // Sticky action
                that.actions = new Sticky({
                    debug: that.config.debug,
                    element: that.$actions,
                    parent: that.$document,
                    top: 20,
                    bottom: 20
                });
                that.actions.enable();

                that.$document.show();

                that.loaded = true;
                that.loadXhr = null;

                $(that).trigger('etalia.detail.loaded');
            });

        return this;
    };

    Detail.prototype.clear = function () {
        this.loaded = false;
        if (this.loadXhr) {
            this.loadXhr.abort();
            this.loadXhr = null;
        }

        if (this.actions) {
            this.actions.disable();
            this.actions = null;
        }

        if (this.$actions) {
            this.$actions = null;
        }

        if (this.paperId) {
            this.paperId = null;
        }

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
