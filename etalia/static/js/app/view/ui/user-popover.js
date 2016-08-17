define([
    'app/app',
    'app/model/ui/user-popover',
    'bootstrap'
], function (App) {

    var _ = App._,
        defaults = {
            $target: null
        };

    App.View.Ui = App.View.Ui || {};

    App.View.Ui.UserPopover = App.Backbone.View.extend({

        $target: null,
        targetPosition: {top: 0, left: 0},
        position: {top: 0, left: 0},

        popover: null,

        attributes: {
            'class': 'ui-user-popover'
        },

        events: {
            'click': 'onClick'
        },

        initialize: function (options) {
            options || (options = {});
            _.defaults(this, defaults);
            _.extend(this, _.pick(options, _.keys(defaults)));

            _.bindAll(this, '_watchPosition', '_updatePosition');

            setTimeout(this._watchPosition, 1000);
            App.$(window).on('resize', this._watchPosition);
            App.$(document).on('scroll', this._watchPosition);
        },

        onClick: function() {

            App.$('.ui-user-popover').popover('hide');

            if (this.popover) {
                this.$el.popover('show');
                return;
            }

            var that = this,
                data = this.model.get('popover');
            this.popover = this.$el.popover({
                title: data.get('title'),
                content: data.get('body'),
                html: true,
                placement: 'auto',
                template: '<div class="popover" role="tooltip">' +
                    '<div class="arrow"></div>' +
                    //'<h3 class="popover-title"></h3>' +
                    '<div class="popover-content"></div>' +
                    '<div class="popover-footer">' +
                        '<button class="btn btn-primary btn-sm">Got it !</button>' +
                    '</div>' +
                '</div>'
            }).data('bs.popover');

            this.$el.popover('show');

            this.popover.$tip.find('.popover-footer .btn').on('click', function() {
                that.model.markDone().then(function() {
                    that.$el.popover('destroy');
                    that.remove();
                });
            });
        },

        render: function () {
            this.$el.html(
                //'<span class="eai eai-exclamation"></span>' +
                '<span class="exclamation">?</span>' +
                '<span class="circle"></span>' +
                '<span class="circle"></span>' +
                '<span class="circle"></span>'
            );

            this._watchPosition();

            return this;
        },

        _watchPosition: function() {
            var targetPosition = this.$target.offset();
            if (targetPosition != this.targetPosition) {
                this.targetPosition = targetPosition;
                this._updatePosition();
            }
        },

        _updatePosition: function() {
            this.position = {top: 0, left: 0};
            this.position.top = this.targetPosition.top + (this.$target.height() / 2) + 12;
            this.position.left = this.targetPosition.left + (this.$target.width() / 2) + 12;
            this.$el.css(this.position);
        },

        _openPopover: function() {
            var popover = this.model.get('popover'),
                options = {
                    title: data.get('title'),
                    content: data.get('content'),
                    html: true
                };

            var windowSize = {
                    width: $(window).width(),
                    height: $(window).height()
                };


            //options.placement;

            this.$el.popover(options);
        }
    });

    App.View.Ui.UserPopover.create = function(options) {
        var view = new App.View.Ui.UserPopover(options);
        $('body').append(view.render().$el);
        return view;
    };

    return App.View.Ui.UserPopover;
});
